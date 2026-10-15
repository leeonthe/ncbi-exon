"""CLI entry point.

Usage:
    python -m ncbi_exon \
        --input "Protein Coding genes (human).xlsx" \
        --output output/results.csv \
        --cache-dir cache \
        [--refresh] [--rps 10] [--concurrency 10] [--taxon human] [--limit N] \
        [--no-protein-coding-filter] [--log output/run.log]

Reads gene symbols from CSV/XLSX, fetches each gene's product report from the
NCBI Datasets v2 API, parses, and writes one CSV row per input gene. Resumable
via the file cache (skip cache with --refresh). Ctrl-C flushes partial output.
"""

from __future__ import annotations

import argparse
import asyncio
import signal
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .async_client import AsyncNCBIClient, FetchResult
from .env import load_dotenv
from .input_reader import read_gene_symbols
from .output_writer import GeneRowCSVWriter
from .parser import GeneRow, merge_paged_payloads, parse_product_report


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m ncbi_exon", description=__doc__)
    p.add_argument("--input", "-i", required=True, help="Input file (.csv or .xlsx)")
    p.add_argument("--output", "-o", required=True, help="Output CSV path")
    p.add_argument("--cache-dir", default="cache", help="API response cache directory")
    p.add_argument("--log", default=None, help="Path to run.log (default: <output>.log)")
    p.add_argument("--taxon", default="human", help="NCBI taxon (default: human)")
    p.add_argument("--rps", type=float, default=None, help="Requests/sec (default: 10 w/ key, 5 w/o)")
    p.add_argument("--concurrency", type=int, default=None, help="Concurrent in-flight requests (default: rps)")
    p.add_argument("--limit", type=int, default=None, help="Process only the first N unique symbols")
    p.add_argument("--refresh", action="store_true", help="Bypass cache, re-fetch every gene")
    p.add_argument("--no-protein-coding-filter", action="store_true",
                   help="Don't filter input rows by type_of_gene == protein-coding")
    p.add_argument("--sheet", default=None, help="XLSX sheet name (default: first sheet)")
    return p


async def _run(args: argparse.Namespace) -> int:
    load_dotenv()

    print(f"[{_now()}] reading input: {args.input}")
    symbols, stats = read_gene_symbols(
        args.input,
        sheet=args.sheet,
        filter_protein_coding=not args.no_protein_coding_filter,
        limit=args.limit,
    )
    print(f"[{_now()}] input scan: total_rows={stats.total_rows}  "
          f"symbol_col={stats.symbol_column!r}  type_col={stats.type_column!r}  "
          f"after_filter={stats.after_filter}  unique={stats.after_dedupe}  "
          f"dropped_empty={stats.dropped_empty}")
    if not symbols:
        print("[error] no symbols to process", file=sys.stderr)
        return 1

    output_path = Path(args.output)
    log_path = Path(args.log) if args.log else output_path.with_suffix(output_path.suffix + ".log")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total = len(symbols)
    counts = {"OK": 0, "NOT_FOUND": 0, "MULTI_MATCH": 0, "ERROR": 0, "is_exon_true": 0}
    cancelled = False

    # Graceful Ctrl-C: set a flag; the producer loop checks it and stops scheduling new work.
    def _on_sigint(*_):
        nonlocal cancelled
        if cancelled:
            return  # second Ctrl-C → let it raise
        cancelled = True
        print(f"\n[{_now()}] SIGINT received — finishing in-flight requests, then flushing...",
              file=sys.stderr)

    signal.signal(signal.SIGINT, _on_sigint)

    print(f"[{_now()}] output:  {output_path}")
    print(f"[{_now()}] log:     {log_path}")
    print(f"[{_now()}] cache:   {args.cache_dir}")

    t_start = time.monotonic()

    with GeneRowCSVWriter(output_path) as writer, log_path.open("w", encoding="utf-8") as log_fh:
        log_fh.write(f"# ncbi-exon run started {_now()}  total={total}\n")
        log_fh.write("# ts | symbol | status | http | pages | latency_ms | nm | nr | is_exon\n")
        log_fh.flush()

        async with AsyncNCBIClient(
            rps=args.rps,
            concurrency=args.concurrency,
            cache_dir=args.cache_dir,
        ) as client:
            print(f"[{_now()}] client: rps={client.rps}  concurrency={client.concurrency}  "
                  f"api_key={'yes' if client.api_key else 'no'}")

            async def _fetch(idx: int, sym: str) -> tuple[int, str, FetchResult | Exception]:
                try:
                    res = await client.fetch_product_report(sym, taxon=args.taxon, refresh=args.refresh)
                    return idx, sym, res
                except Exception as e:  # noqa: BLE001 — capture any to keep the pipeline alive
                    return idx, sym, e

            # In-order writer: tasks complete out-of-order, but we buffer results in a
            # dict keyed by input index and flush a contiguous prefix to the CSV so the
            # output is byte-identical to the input order regardless of latency variance.
            buffer: dict[int, tuple[str, object, GeneRow]] = {}
            next_write_idx = 0
            pending: set[asyncio.Task] = set()
            sym_iter = enumerate(symbols)
            max_inflight = client.concurrency * 4
            done_count = 0

            def _schedule_more() -> None:
                while not cancelled and len(pending) < max_inflight:
                    try:
                        i, s = next(sym_iter)
                    except StopIteration:
                        return
                    pending.add(asyncio.create_task(_fetch(i, s)))

            _schedule_more()

            while pending:
                done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
                for task in done:
                    idx, sym, res = task.result()
                    row = _result_to_row(sym, res)
                    buffer[idx] = (sym, res, row)

                # Flush contiguous prefix.
                while next_write_idx in buffer:
                    sym, res, row = buffer.pop(next_write_idx)
                    writer.write(row)
                    counts[row.status] = counts.get(row.status, 0) + 1
                    if row.is_exon:
                        counts["is_exon_true"] += 1
                    done_count += 1
                    log_fh.write(
                        f"{_now()} | {sym} | {row.status} | "
                        f"{getattr(res, 'status_code', 0) if not isinstance(res, Exception) else 0} | "
                        f"{getattr(res, 'pages', 0) if not isinstance(res, Exception) else 0} | "
                        f"{getattr(res, 'latency_ms', 0) if not isinstance(res, Exception) else 0} | "
                        f"{row.nm_count} | {row.nr_count} | {row.is_exon}\n"
                    )
                    log_fh.flush()
                    next_write_idx += 1
                    if done_count % 25 == 0 or done_count == total:
                        elapsed = time.monotonic() - t_start
                        rate = done_count / elapsed if elapsed > 0 else 0
                        eta = (total - done_count) / rate if rate > 0 else 0
                        print(f"  [{done_count:>5}/{total}]  OK={counts['OK']}  "
                              f"NOT_FOUND={counts['NOT_FOUND']}  MULTI={counts['MULTI_MATCH']}  "
                              f"ERROR={counts['ERROR']}  is_exon=True×{counts['is_exon_true']}  "
                              f"{rate:.1f} req/s  ETA {eta/60:.1f} min")

                _schedule_more()

    elapsed = time.monotonic() - t_start
    print()
    print(f"[{_now()}] DONE  in {elapsed:.1f}s  ({done_count}/{total} rows)")
    print(f"  OK={counts['OK']}  NOT_FOUND={counts['NOT_FOUND']}  "
          f"MULTI_MATCH={counts['MULTI_MATCH']}  ERROR={counts['ERROR']}")
    print(f"  is_exon=True: {counts['is_exon_true']}")
    print(f"  output: {output_path}")
    print(f"  log:    {log_path}")
    if cancelled:
        print("  [partial run — interrupted by user]")
        return 130
    return 0


def _result_to_row(symbol: str, res) -> GeneRow:
    """Convert a FetchResult (or Exception) to a GeneRow."""
    if isinstance(res, Exception):
        return parse_product_report(symbol, payload={}, error_detail=f"{type(res).__name__}: {res}")
    if res.error_detail:
        return parse_product_report(symbol, payload={}, error_detail=res.error_detail)
    merged = merge_paged_payloads(res.payloads) if res.payloads else {}
    row = parse_product_report(symbol, merged, http_status=res.status_code)
    return row


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    try:
        return asyncio.run(_run(args))
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())

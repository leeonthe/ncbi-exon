"""Async NCBI Datasets v2 client (httpx).

Used by the bulk CLI runner. Honors the global rate limit (10 rps with API key,
5 rps without), runs requests concurrently via a semaphore, retries 429/5xx with
exponential backoff, and follows `next_page_token` so multi-page responses are
fully fetched before the caller sees them.

The synchronous `client.NCBIClient` is kept for the probe script and small ad-hoc
tests; this module is for the 20k-gene run.
"""

from __future__ import annotations

import asyncio
import json
import os
import time
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import httpx

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_TIMEOUT = 30.0
DEFAULT_MAX_RETRIES = 4
RETRYABLE_STATUS = {429, 500, 502, 503, 504}
MAX_PAGES = 20  # defensive cap — a single gene should never produce this many pages


class NCBIError(Exception):
    def __init__(self, message: str, status: int | None = None):
        super().__init__(message)
        self.status = status


@dataclass
class FetchResult:
    symbol: str
    status_code: int
    payloads: list[dict[str, Any]]   # one entry per page; len==1 in the common case
    latency_ms: int
    from_cache: bool = False
    pages: int = 1
    error_detail: str = ""


class _RateLimiter:
    """Simple async leaky-bucket: ensures gaps of >= 1/rps between request starts.

    With a semaphore of size N >= rps, this is sufficient to keep us right at the
    NCBI rate limit without bursting over it.
    """

    def __init__(self, rps: float):
        self._min_interval = 1.0 / rps if rps > 0 else 0.0
        self._next_at = 0.0
        self._lock = asyncio.Lock()

    async def gate(self) -> None:
        async with self._lock:
            now = time.monotonic()
            wait = self._next_at - now
            if wait > 0:
                await asyncio.sleep(wait)
                now = time.monotonic()
            self._next_at = max(self._next_at, now) + self._min_interval


class AsyncNCBIClient:
    """Async client. Use within an `async with` block."""

    def __init__(
        self,
        *,
        api_key: str | None = None,
        rps: float | None = None,
        concurrency: int | None = None,
        timeout: float = DEFAULT_TIMEOUT,
        max_retries: int = DEFAULT_MAX_RETRIES,
        cache_dir: str | Path | None = None,
    ):
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.rps = rps if rps is not None else (10.0 if self.api_key else 5.0)
        self.concurrency = concurrency if concurrency is not None else max(int(self.rps), 1)
        self.timeout = timeout
        self.max_retries = max_retries
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._limiter = _RateLimiter(self.rps)
        self._sem = asyncio.Semaphore(self.concurrency)
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self) -> "AsyncNCBIClient":
        headers = {"Accept": "application/json", "User-Agent": "ncbi-exon/0.1"}
        if self.api_key:
            headers["api-key"] = self.api_key
        self._client = httpx.AsyncClient(timeout=self.timeout, headers=headers)
        return self

    async def __aexit__(self, *exc) -> None:
        if self._client:
            await self._client.aclose()
            self._client = None

    async def fetch_product_report(
        self,
        symbol: str,
        taxon: str = "human",
        *,
        refresh: bool = False,
    ) -> FetchResult:
        """Fetch the full (paginated) product report for one symbol."""
        symbol = symbol.strip()
        if not symbol:
            raise ValueError("symbol must be non-empty")

        cache_path = self._cache_path(symbol, taxon)
        if cache_path and not refresh and cache_path.exists():
            try:
                payload = json.loads(cache_path.read_text())
                return FetchResult(symbol=symbol, status_code=200, payloads=[payload],
                                   latency_ms=0, from_cache=True, pages=1)
            except json.JSONDecodeError:
                pass  # corrupt cache, refetch

        first_url = (
            f"{API_BASE}/gene/symbol/{urllib.parse.quote(symbol, safe='')}"
            f"/taxon/{urllib.parse.quote(taxon, safe='')}/product_report"
        )

        t0 = time.monotonic()
        payloads: list[dict[str, Any]] = []
        url: str | None = first_url
        status = 0
        page_count = 0
        err_detail = ""

        async with self._sem:
            while url and page_count < MAX_PAGES:
                page_count += 1
                status, body, err = await self._get(url)
                if err:
                    err_detail = err
                    break
                if status != 200:
                    # 404 or other 4xx — payload may carry an error envelope; we still want it.
                    if body:
                        try:
                            payloads.append(json.loads(body))
                        except json.JSONDecodeError:
                            pass
                    break
                try:
                    pl = json.loads(body)
                except json.JSONDecodeError as e:
                    err_detail = f"invalid JSON: {e}"
                    break
                payloads.append(pl)
                token = pl.get("next_page_token")
                if not token:
                    break
                url = (
                    f"{first_url}?page_token={urllib.parse.quote(token, safe='')}"
                )

        latency_ms = int((time.monotonic() - t0) * 1000)

        # Cache only on a single-page 200 (multi-page caching is unnecessary in practice).
        if status == 200 and not err_detail and cache_path and len(payloads) == 1:
            try:
                cache_path.write_text(json.dumps(payloads[0]))
            except OSError:
                pass

        return FetchResult(
            symbol=symbol,
            status_code=status,
            payloads=payloads,
            latency_ms=latency_ms,
            from_cache=False,
            pages=page_count,
            error_detail=err_detail,
        )

    # ------------------------------------------------------------ internals

    def _cache_path(self, symbol: str, taxon: str) -> Path | None:
        if not self.cache_dir:
            return None
        safe = symbol.upper().replace("/", "_")
        return self.cache_dir / f"{safe}__{taxon}.json"

    async def _get(self, url: str) -> tuple[int, bytes, str]:
        assert self._client is not None
        backoff = 1.0
        last_err = ""
        for attempt in range(self.max_retries + 1):
            await self._limiter.gate()
            try:
                resp = await self._client.get(url)
                if resp.status_code in RETRYABLE_STATUS and attempt < self.max_retries:
                    retry_after = self._parse_retry_after(resp.headers.get("Retry-After"))
                    await asyncio.sleep(retry_after if retry_after is not None else backoff)
                    backoff = min(backoff * 2, 16.0)
                    continue
                return resp.status_code, resp.content, ""
            except (httpx.TimeoutException, httpx.TransportError) as e:
                last_err = f"{type(e).__name__}: {e}"
                if attempt < self.max_retries:
                    await asyncio.sleep(backoff)
                    backoff = min(backoff * 2, 16.0)
                    continue
                return 0, b"", last_err
        return 0, b"", last_err or "exhausted retries"

    @staticmethod
    def _parse_retry_after(value: str | None) -> float | None:
        if not value:
            return None
        try:
            return float(value)
        except ValueError:
            return None

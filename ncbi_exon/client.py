"""HTTP client for the NCBI Datasets v2 REST API.

Validation-phase implementation: stdlib only (urllib), synchronous, polite pacing
between requests. Swap to async httpx when scaling past a few hundred genes.

Endpoint used:
    GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/{taxon}/product_report

See README.md §2 for the API decision and §6 for the runtime requirements.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_TIMEOUT = 30.0
DEFAULT_MAX_RETRIES = 4
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class NCBIError(Exception):
    """Raised for unrecoverable NCBI API errors (after retries are exhausted)."""

    def __init__(self, message: str, status: int | None = None):
        super().__init__(message)
        self.status = status


@dataclass
class FetchResult:
    """Result of a single product_report fetch."""

    symbol: str
    status_code: int           # HTTP status of the final attempt (200 on success, 404 etc.)
    payload: dict[str, Any]    # parsed JSON body (empty dict on non-200)
    latency_ms: int            # total wall time including retries
    from_cache: bool = False


class NCBIClient:
    """Synchronous client for NCBI Datasets v2 gene product reports.

    Reads API key from the NCBI_API_KEY env var (or pass api_key=). Paces requests
    at `rps` per second using a simple sleep-based gate.
    """

    def __init__(
        self,
        *,
        api_key: str | None = None,
        rps: float | None = None,
        timeout: float = DEFAULT_TIMEOUT,
        max_retries: int = DEFAULT_MAX_RETRIES,
        cache_dir: str | Path | None = None,
    ):
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        # 10 rps with key, 5 rps without — per NCBI docs
        self.rps = rps if rps is not None else (10.0 if self.api_key else 5.0)
        self.timeout = timeout
        self.max_retries = max_retries
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._min_interval = 1.0 / self.rps if self.rps > 0 else 0.0
        self._last_request_at = 0.0

    # ------------------------------------------------------------------ public

    def fetch_product_report(
        self,
        symbol: str,
        taxon: str = "human",
        *,
        refresh: bool = False,
    ) -> FetchResult:
        """Fetch the gene product report for a single symbol+taxon.

        Returns a FetchResult. Does NOT raise on 404 / empty results — those are
        legitimate outcomes recorded by the caller. Raises NCBIError only when
        retries are exhausted on a transient failure.
        """
        symbol = symbol.strip()
        if not symbol:
            raise ValueError("symbol must be non-empty")

        cache_path = self._cache_path(symbol, taxon)
        if cache_path and not refresh and cache_path.exists():
            payload = json.loads(cache_path.read_text())
            return FetchResult(symbol=symbol, status_code=200, payload=payload, latency_ms=0, from_cache=True)

        url = (
            f"{API_BASE}/gene/symbol/{urllib.parse.quote(symbol, safe='')}"
            f"/taxon/{urllib.parse.quote(taxon, safe='')}/product_report"
        )
        headers = {"Accept": "application/json", "User-Agent": "ncbi-exon/0.1"}
        if self.api_key:
            headers["api-key"] = self.api_key

        t0 = time.monotonic()
        status, body = self._request_with_retries(url, headers)
        latency_ms = int((time.monotonic() - t0) * 1000)

        payload: dict[str, Any] = {}
        if body:
            try:
                payload = json.loads(body)
            except json.JSONDecodeError as e:
                raise NCBIError(f"Invalid JSON from {url}: {e}") from e

        if status == 200 and cache_path:
            cache_path.write_text(json.dumps(payload, indent=2))

        return FetchResult(symbol=symbol, status_code=status, payload=payload, latency_ms=latency_ms)

    # ----------------------------------------------------------------- internal

    def _cache_path(self, symbol: str, taxon: str) -> Path | None:
        if not self.cache_dir:
            return None
        safe = symbol.upper().replace("/", "_")
        return self.cache_dir / f"{safe}__{taxon}.json"

    def _gate(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_at = time.monotonic()

    def _request_with_retries(self, url: str, headers: dict[str, str]) -> tuple[int, bytes]:
        backoff = 1.0
        last_exc: Exception | None = None
        for attempt in range(self.max_retries + 1):
            self._gate()
            try:
                req = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    return resp.status, resp.read()
            except urllib.error.HTTPError as e:
                # Read body for diagnostics (NCBI returns JSON error envelopes).
                body = b""
                try:
                    body = e.read()
                except Exception:
                    pass
                if e.code in RETRYABLE_STATUS and attempt < self.max_retries:
                    retry_after = self._parse_retry_after(e.headers.get("Retry-After"))
                    time.sleep(retry_after if retry_after is not None else backoff)
                    backoff = min(backoff * 2, 16.0)
                    continue
                # 4xx (other than 429) is terminal: surface status + body to caller.
                return e.code, body
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                last_exc = e
                if attempt < self.max_retries:
                    time.sleep(backoff)
                    backoff = min(backoff * 2, 16.0)
                    continue
                raise NCBIError(f"Transport error after {self.max_retries} retries: {e}") from e
        # Unreachable, but keeps type checkers happy.
        raise NCBIError(f"Exhausted retries: {last_exc}")

    @staticmethod
    def _parse_retry_after(value: str | None) -> float | None:
        if not value:
            return None
        try:
            return float(value)
        except ValueError:
            return None

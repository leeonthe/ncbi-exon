"""Streaming CSV writer for GeneRow records.

Writes one row per gene with a fixed column order matching README.md §5.1
(minus `synonyms` and `ensembl_gene_id` per Option A decision 2026-06-02).

Use as a context manager:
    with GeneRowCSVWriter(path) as w:
        for row in rows:
            w.write(row)

The writer flushes after every write, so partial output survives Ctrl-C.
"""

from __future__ import annotations

import csv
from dataclasses import fields
from pathlib import Path
from types import TracebackType
from typing import IO

from .parser import GeneRow

# Fixed column order — must match README §5.1 (minus synonyms/ensembl).
COLUMNS = [f.name for f in fields(GeneRow)]


class GeneRowCSVWriter:
    """Streaming CSV writer. Flushes after every row to survive Ctrl-C."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._fh: IO[str] | None = None
        self._writer: csv.DictWriter | None = None

    def __enter__(self) -> "GeneRowCSVWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self.path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(self._fh, fieldnames=COLUMNS)
        self._writer.writeheader()
        self._fh.flush()
        return self

    def write(self, row: GeneRow) -> None:
        assert self._writer is not None and self._fh is not None
        d = row.as_dict()
        # csv.DictWriter writes booleans as 'True'/'False' — keep that, it round-trips fine.
        self._writer.writerow(d)
        self._fh.flush()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._fh:
            self._fh.flush()
            self._fh.close()
            self._fh = None
            self._writer = None

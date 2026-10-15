"""NCBI Exon — fetch RefSeq transcript counts for human genes via NCBI Datasets v2."""

from .client import NCBIClient, NCBIError
from .parser import GeneRow, merge_paged_payloads, parse_product_report

__all__ = [
    "NCBIClient",
    "NCBIError",
    "parse_product_report",
    "merge_paged_payloads",
    "GeneRow",
]

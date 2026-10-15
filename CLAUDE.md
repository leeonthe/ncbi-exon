# Project rules for Claude Code

**MANDATORY: Before writing or modifying ANY code in this repository, read `README.md` in full.**

`README.md` is the source of truth for:
- The project goal (NCBI gene → transcript annotation for ~20k human genes)
- The chosen API (NCBI Datasets v2 REST — endpoint, headers, response shape)
- Counting rules (NM/NR + v2: XM/XR; `is_exon = (nm_count == 1)` strict equality)
- Input/output schema and column definitions (v2 schema: 32 columns — 29 base + 3 Definition-#2 containment columns)
- Rate limits, retries, caching, resumability requirements
- §1.5 — current implementation status (what is done and what is left)
- §8 — verification checklist

Companion docs (load when their topic is in scope):
- `Definition_of_non_unique_transcript.md` — user-supplied definition of the overlap rule
- `Approach_non_unique_transcripts.md` — logic spec for v2 overlap analysis, with worked examples (ACTB / MYC / GAPDH) computed from real API data
- `Visualization_plan.md` — per-figure rationale, statistical standards, thesis narrative, Phase-2 figure proposals
- `How_to_generate_figures.md` — paste-into-Jupyter snippets per figure with overlap-free layouts; Plotly Sankey alternative for the flow diagram

If a user request seems to conflict with `README.md`, surface the conflict and ask before changing code. Update `README.md` first when requirements change, then update the code to match — never the other way around.

Other constraints:
- `NCBI_API_KEY` is read from env. Never print, log, or commit it.
- Do not introduce XM/XR counting unless `README.md` is updated to include it (currently out of scope).
- Cache every API response to `cache/<SYMBOL>.json` so reruns don't re-hit the API.

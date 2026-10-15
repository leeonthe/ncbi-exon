# Approach — Non-Unique Transcript Detection (v2 output)

This is the **logic** spec for the next rebuild of `output.csv`. It does not describe code — it describes what we will compute, what data each value is based on, and walks through three real genes using the actual NCBI API response so you can confirm the rule matches your intent.

Companion to `Definition_of_non_unique_transcript.md` (the user-supplied definition).

---

## 1. What we want to know about each gene

For every gene in the input we already record `is_exon`, `nm_count`, `nr_count`, and the supporting accessions. **v2 adds two new things:**

1. **Predicted transcript counts** — how many `XM_` (predicted protein-coding) and `XR_` (predicted non-coding) RefSeq transcripts exist, alongside the curated `NM_` / `NR_`.
2. **Non-unique transcript flag** — does this gene contain at least one transcript whose every exon overlaps an exon of some other transcript? (Per the definition in `Definition_of_non_unique_transcript.md`.)

---

## 2. What values we'll store, and where each comes from

All values come from a **single API call per gene**:

```
GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/human/product_report
```

That response gives us, per transcript: the `accession_version` (e.g. `NM_002467.6`), the genomic assembly it maps to (e.g. GRCh38), and the list of exons on that assembly with start/end coordinates. Everything below is derived from those fields — no extra API calls.

### Existing columns (unchanged)

| Column | Meaning | Source |
|---|---|---|
| `nm_count`, `nr_count` | counts of `NM_…` / `NR_…` accessions | parse `transcripts[].accession_version` |
| `is_exon` | `True` iff `nm_count == 1` | derived |
| `mane_select_accession` | the MANE Select transcript, if any | `transcripts[].select_category == "MANE_SELECT"` |
| ... | (all other existing fields) | |

### New columns

| Column | Type | Meaning | Source |
|---|---|---|---|
| `xm_count` | int | predicted protein-coding transcripts | accessions starting with `XM_` |
| `xr_count` | int | predicted non-coding transcripts | accessions starting with `XR_` |
| `xm_accessions` | str | `;`-joined `XM_` list | |
| `xr_accessions` | str | `;`-joined `XR_` list | |
| `has_non_unique_transcript` | bool | gene-level overlap verdict | derived from exon coordinates (see §3) |
| `non_unique_transcript_count` | int | how many transcripts of this gene are non-unique | |
| `non_unique_transcript_accessions` | str | `;`-joined list of those transcripts | |
| `transcripts_compared_for_overlap` | int | how many transcripts had usable GRCh38 exon data and were included in the analysis | |
| `transcripts_excluded_from_overlap` | int | how many were skipped (no GRCh38 location, no exons, etc.) | |

---

## 3. The non-unique transcript rule (in plain English)

This is the rule from your definition file, restated to be unambiguous:

> A transcript **T** is **non-unique** if there exists *some other* transcript **T′** (of the same gene) such that **every** exon of **T** overlaps with **at least one** exon of **T′**.
>
> A gene **has a non-unique transcript** if **at least one** of its transcripts is non-unique.

Three things to notice:

1. The check is **directional**: "T is non-unique vs T′" is not the same as "T′ is non-unique vs T". (Your Gene Y example shows this: T1 is non-unique vs T3, but T3 has a unique exon vs T1.)
2. We compare **pairwise** — one transcript vs one other transcript at a time. We do **not** take the union of all other transcripts' exons.
3. We declare the gene "has non-unique transcript" the moment **any** transcript turns out to be non-unique vs **any** other. We don't require all of them to be.

### What "overlap" means for two exons

Two exons overlap if their genomic ranges share **at least one base pair** on the same chromosome:

```
Exon A: 100–200       Exon B: 150–300       → overlap (share 150–200)
Exon A: 100–200       Exon B: 201–300       → no overlap (touching but not sharing)
Exon A: 100–200       Exon B: 50–250        → overlap (A is fully inside B)
```

Standard genomic intervals, 1-based inclusive (which is what NCBI gives us).

### What assembly do we use?

Each NCBI transcript lists exons against multiple assemblies (e.g. GRCh38 *and* T2T-CHM13). **We will only use GRCh38** (the reference standard, accession prefix `NC_0000…`) and drop the alternates. If we mixed assemblies, the coordinates wouldn't line up and overlap would be meaningless.

### Which transcripts go into the comparison?

**`NM_` (curated protein-coding mRNA) only.** Updated 2026-06-02 — predicted transcripts (`XM_`, `XR_`) are treated as "not yet confirmed" and curated non-coding RNAs (`NR_`) are out of scope for the coding-isoform redundancy question. The goal of this analysis is whether the *protein-coding* annotation contains redundant isoforms; that is a question about curated mRNAs.

Concretely, this means:

- `XM_` and `XR_` are reported in their dedicated count columns (`xm_count`, `xr_count`) but never participate in the overlap test.
- `NR_` is reported in `nr_count` but does not participate either.
- An NM transcript with no usable GRCh38 location or no exons is dropped and counted in `transcripts_excluded_from_overlap`.
- `transcripts_compared_for_overlap` always equals "NMs that had usable GRCh38 data and were tested pairwise."

This change makes `has_non_unique_transcript` a statement about the curated protein-coding annotation specifically, which is what the project's downstream use (isoform-specific assays, primer design) actually cares about.

---

## 3.1 Definition #2 — the containment ("2-rule") variant

The rule in §3 above is **Definition #1** (the "3-rule" definition): two exons overlap if they share **at least one base pair** (partial overlap counts). The three illustrative cases are:

```
Exon A: 100–200   Exon B: 150–300   → overlap (share 150–200)   ← partial overlap counts
Exon A: 100–200   Exon B: 201–300   → no overlap (touching but not sharing)
Exon A: 100–200   Exon B: 50–250    → overlap (A is fully inside B)
```

**Definition #2** (the "2-rule" definition) keeps only the **last two** cases and drops the partial-overlap case. Under Definition #2, two exons only count as overlapping when **one is fully contained inside the other** — a partial overlap is treated as *no overlap*:

```
Exon A: 100–200   Exon B: 201–300   → no overlap (touching but not sharing)
Exon A: 100–200   Exon B: 50–250    → overlap (A is fully inside B)
Exon A: 100–200   Exon B: 150–300   → NO overlap (only partial — does not count under #2)
```

**Containment is directional.** Restating the non-unique rule for Definition #2:

> A transcript **T** is **non-unique (containment)** if there exists *some other* transcript **T′** (of the same gene) such that **every** exon of **T** lies **fully inside** at least one exon of **T′** (`exon.begin ≥ T′_exon.begin` AND `exon.end ≤ T′_exon.end`).

So "T is absorbed into T′" — every exon of T is a sub-interval of one of T′'s exons. The reverse direction (an exon of T containing an exon of T′) does **not** satisfy the candidate's requirement.

Both definitions share everything else: the NM-only candidate set, GRCh38-only coordinates, pairwise one-vs-one comparison (never the union), and the "any transcript flips the gene" rule. Only the per-pair predicate changes.

**Relationship between the two:** containment implies overlap, so **Definition #2 is strictly stronger** — any transcript flagged non-unique under #2 is also flagged under #1, but not vice-versa. The #2 flagged set is always a subset of the #1 flagged set.

### Output columns

Both definitions are computed on every run and stored side by side. Definition #1 keeps the original column names; Definition #2 adds three columns with a `_containment` suffix:

| Column | Definition | Meaning |
|---|---|---|
| `has_non_unique_transcript` | #1 partial | gene has ≥1 non-unique NM (partial-overlap rule) |
| `non_unique_transcript_count` | #1 partial | how many NMs are non-unique under #1 |
| `non_unique_transcript_accessions` | #1 partial | `;`-joined list |
| `has_non_unique_transcript_containment` | #2 containment | gene has ≥1 non-unique NM (containment rule) |
| `non_unique_transcript_count_containment` | #2 containment | how many NMs are non-unique under #2 |
| `non_unique_transcript_accessions_containment` | #2 containment | `;`-joined list |
| `transcripts_compared_for_overlap` | shared | NM candidate count (same for both definitions) |
| `transcripts_excluded_from_overlap` | shared | NMs dropped for missing GRCh38 data (same for both) |

### Worked-example deltas (real cached data)

| Gene | #1 count | #2 count | Note |
|---|---|---|---|
| ACTB | 0 | 0 | single NM — trivially clean under both |
| MYC | 2 | **0** | both NMs partially overlap, but neither is *fully contained* (T1's exon 1 extends ~800 bp past T2's; T2's exon 2 starts 3 bp before T1's) — so the stricter rule flags nothing |
| GAPDH | 4 | 3 | NM_001289746.2 drops out under containment |
| BRCA2 | 5 | 4 | NM_000059.4 (the MANE Select) drops out under containment |

MYC is the cleanest illustration of why the rules differ: two transcripts can be "the same gene body with a slightly longer first exon," which is redundant under partial overlap but distinguishable (not contained) under the stricter rule.

---

## 4. Worked example #1 — ACTB (simple, "no non-unique")

**Input:** ACTB (gene_id 60), chromosome 7. The API returns:

| Transcript | GRCh38 location | Exons (begin–end) |
|---|---|---|
| NM_001101.5 | NC_000007.14 | 5527148–5527891 · 5528004–5528185 · 5528281–5528719 · 5529161–5529400 · 5529535–5529663 · 5530524–5530601 |

Only **one** transcript exists with usable exon data. There is nothing to compare against.

**Result:**
- `xm_count = 0`, `xr_count = 0`, `nm_count = 1`, `nr_count = 0`
- `transcripts_compared_for_overlap = 1`
- `non_unique_transcript_count = 0`
- **`has_non_unique_transcript = False`** (no other transcript exists, so by definition no transcript can be non-unique).

---

## 5. Worked example #2 — MYC (positive, "has non-unique")

**Input:** MYC (gene_id 4609), chromosome 8. The API returns 2 NM transcripts, both on GRCh38:

| Transcript | Exon 1 | Exon 2 | Exon 3 |
|---|---|---|---|
| **NM_001354870.1** (T1) | 127,735,434–127,736,623 | 127,738,251–127,739,019 | 127,740,396–127,742,951 |
| **NM_002467.6** (T2, MANE Select) | 127,736,231–127,736,623 | 127,738,248–127,739,019 | 127,740,396–127,742,951 |

These transcripts are nearly identical — the only real difference is that T1's first exon starts ~800 bp earlier than T2's first exon. Exons 2 and 3 are essentially the same.

### Pairwise check: is T1 non-unique vs T2?

We ask: does **every** exon of T1 overlap **at least one** exon of T2?

| T1 exon | overlaps any T2 exon? |
|---|---|
| Exon 1 (127,735,434–127,736,623) | ✅ overlaps T2 Exon 1 (127,736,231–127,736,623) — they share positions 127,736,231–127,736,623 |
| Exon 2 (127,738,251–127,739,019) | ✅ overlaps T2 Exon 2 |
| Exon 3 (127,740,396–127,742,951) | ✅ overlaps T2 Exon 3 |

All three exons of T1 overlap something in T2 → **T1 is non-unique vs T2.**

### Pairwise check: is T2 non-unique vs T1?

| T2 exon | overlaps any T1 exon? |
|---|---|
| Exon 1 (127,736,231–127,736,623) | ✅ overlaps T1 Exon 1 |
| Exon 2 (127,738,248–127,739,019) | ✅ overlaps T1 Exon 2 |
| Exon 3 (127,740,396–127,742,951) | ✅ overlaps T1 Exon 3 |

Also yes → **T2 is non-unique vs T1.**

**Result:**
- `xm_count = 0`, `xr_count = 0`, `nm_count = 2`, `nr_count = 0`
- `transcripts_compared_for_overlap = 2`
- `non_unique_transcript_count = 2`
- `non_unique_transcript_accessions = "NM_001354870.1;NM_002467.6"`
- **`has_non_unique_transcript = True`**

Note that *both* transcripts are non-unique here, but the gene-level flag would have flipped to `True` as soon as the first one was.

---

## 6. Worked example #3 — GAPDH (complex, NM-only)

**Input:** GAPDH (gene_id 2597), chromosome 12. The API returns 6 transcripts (5 NM + 1 NR), all on GRCh38. **Under the NM-only rule, the single NR transcript is excluded from the comparison and only the 5 NMs are tested pairwise.**

| # | Transcript | First exon | Total exons | In comparison? |
|---|---|---|---|---|
| T1 | NM_001256799.3 | 6,535,264–6,535,469 | 8 | ✅ |
| T2 | NM_001289745.3 | 6,534,517–6,534,661 | 9 | ✅ |
| T3 | NM_001289746.2 | 6,534,517–6,534,861 | 8 | ✅ |
| T4 | NM_001357943.2 | 6,534,517–6,534,569 | 9 | ✅ |
| T5 | NM_002046.7 (MANE) | 6,534,517–6,534,569 | 9 | ✅ |
| — | NR_152150.2 | 6,534,517–6,534,569 | 7 | ❌ (NR, skipped) |

The transcripts mostly share the body of the gene but differ in the 5′ end and one internal exon. We check each NM transcript against every other NM transcript.

After running the pairwise rule:

| Transcript | Result | Why |
|---|---|---|
| T1 (NM_001256799.3) | **unique** | T1's first exon starts at 6,535,264, but every other NM transcript has a 5′ exon farther upstream (6,534,517–6,534,861) that T1 doesn't overlap. T1 has no partner T′ where *all* of T1's exons land inside T′'s exon set. |
| T2 (NM_001289745.3) | **non-unique** | Every exon of T2 overlaps at least one exon in T3 (which has the same 5′ end as a single merged exon plus the rest of the body). |
| T3 (NM_001289746.2) | **non-unique** | Every exon of T3 overlaps at least one exon in T5. |
| T4 (NM_001357943.2) | **non-unique** | T4 differs from T5 by exactly one internal exon coordinate (6,536,738–6,536,790 vs 6,536,684–6,536,790), and that range overlaps T5's exon. Every other exon is identical. |
| T5 (NM_002046.7) | **non-unique** | Symmetric to T4 — every T5 exon overlaps a T4 exon. |

**Result:**
- `xm_count = 0`, `xr_count = 0`, `nm_count = 5`, `nr_count = 1`
- `transcripts_compared_for_overlap = 5` *(NMs only; the NR is not counted)*
- `non_unique_transcript_count = 4`
- `non_unique_transcript_accessions = "NM_001289745.3;NM_001289746.2;NM_001357943.2;NM_002046.7"`
- **`has_non_unique_transcript = True`**

This is a representative real-world case: most multi-NM genes contain at least one "almost-duplicate" curated mRNA pair.

---

## 7. What this flag actually buys you

| If `has_non_unique_transcript = False` | Then |
|---|---|
| Every annotated transcript has at least one exon that you can use to *uniquely* identify it from any other transcript of the same gene. | Good news for primer / probe design, isoform-specific PCR, RNA-seq quantification: every transcript is distinguishable on at least one exon. |

| If `has_non_unique_transcript = True` | Then |
|---|---|
| At least one transcript is genomically a "subset" of another (every exon falls within an exon of another transcript). | Distinguishing those isoforms requires more than just genomic exon presence — you'd need junction-level reads or sequence-level differences (e.g. exon-boundary positions). |

The flag does **not** tell you:

- *Which pair* causes the non-uniqueness (the accession list narrows it but not the partner). We could add that later if needed.
- Whether transcripts differ in CDS vs UTR; the flag is purely exon-coordinate based.
- Whether the difference is biologically meaningful (an exon-boundary shift of 2 bp still counts as overlap by this rule).

---

## 8. Edge cases — how they resolve

| Situation | Behavior |
|---|---|
| Gene has only 1 NM transcript (e.g. ACTB) | `transcripts_compared = 1`, `non_unique_count = 0`, `has_non_unique = False` (trivially — no other NM to compare against). |
| Gene has 0 NM transcripts (e.g. ncRNA-only loci, predicted-only genes) | `transcripts_compared = 0`, `has_non_unique = False`. The XM / NR / XR are reported in their count columns but do not participate. |
| Gene has 0 transcripts at all (NOT_FOUND) | All overlap columns empty; `has_non_unique = False`. |
| An NM transcript has no GRCh38 `genomic_locations` entry | Excluded from the comparison; counted in `transcripts_excluded_from_overlap`. |
| An NM transcript has 0 exons in its GRCh38 entry | Same as above — excluded. |
| All NMs excluded for the above reasons | `transcripts_compared = 0`, `has_non_unique = False`. |
| Transcripts on different strands of the same locus | Compared anyway — overlap is purely coordinate-based. (In practice, all NMs of a single Gene ID are on the same strand, so this rarely matters.) |
| An exon of length 1 (begin == end) | Treated as a 1-bp interval; overlap rule applies normally. |

---

## 9. What we will not do in v2

- We will not switch off the `xm_count` / `xr_count` columns based on the user's earlier "out of scope" decision for `is_exon`. The `is_exon` rule itself stays `nm_count == 1`. The XM/XR columns are pure reporting; they don't change any other value.
- We will not call any additional NCBI endpoint. All v2 work reuses the cached `product_report` payloads.
- We will not invalidate the existing cache. The previous 20,593-row CSV remains valid for the v1 columns; v2 is purely additive.

---

## 10. Test plan

1. Re-run on the first **50 protein-coding genes** from the input file (cache will serve all of them; no API calls needed).
2. Sanity-check expected values (NM-only rule):
   - ACTB → `has_non_unique = False` (single NM)
   - MYC → `has_non_unique = True`, `non_unique_count = 2`
   - GAPDH → `has_non_unique = True`, `non_unique_count = 4` (NR_152150.2 excluded)
   - BRCA2 → `has_non_unique = True`, `non_unique_count = 5` of 6 NMs (NM_001432077.1 is unique — different 5′ exon)
3. Confirm the new columns appear and are populated for `OK` rows, and are empty for `NOT_FOUND` / `MULTI_MATCH` rows.
4. Spot-check a single-transcript gene from the 50 (likely several exist in the first 50) to confirm it lands as `False`.

If those four checks pass, regenerate the full 20,593-row CSV from cache (a few seconds, no API).

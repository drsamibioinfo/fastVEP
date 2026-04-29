# ACMG-AMP Variant Classification in fastVEP: Methods and Benchmark

## Overview

fastVEP implements the 28 ACMG-AMP evidence criteria from Richards et al.
2015 plus published ClinGen Sequence Variant Interpretation (SVI) Working
Group refinements, producing a 5-tier classification: Pathogenic (P),
Likely Pathogenic (LP), Uncertain Significance (VUS), Likely Benign (LB),
Benign (B).

This document reflects the state after the SVI alignment series (PR1–PR10).
Each criterion section includes the ClinGen reference and any deviations
from a strict reading of the SVI text. A per-criterion "Limitations"
column flags criteria that fall back to legacy or conservative behavior
until additional pipeline data is wired in.

## Implementation

### Criteria Coverage

Of the 28 ACMG-AMP criteria, 18 are fully automatable from variant-level
data and are implemented in fastVEP. The classifier records the source
that drove each call in `details.pp3_source` / `details.ps1_path` /
`details.inheritance_basis` etc., so every classification is auditable.

#### Pathogenic Criteria (11 automated)

| Criterion | Strength | Description | Data Source / Notes |
|-----------|----------|-------------|---------------------|
| PVS1 | VS / Strong / Moderate / Supporting (Abou Tayoun 2018 decision tree) | Null variant (nonsense, frameshift, canonical splice, start-loss, whole-gene deletion) in LOF-intolerant gene | Consequence + gnomAD constraints + transcript NMD prediction + critical-region check + alt-start distance |
| PS1 | Strong | Same amino acid change as known pathogenic missense, **or** same RNA outcome for canonical splice (Walker 2023) | ClinVar protein-position index + same-position pathogenic splice catalog |
| PS2 | Strong | Confirmed de novo (full trio) | VCF genotype (proband + both parents) + DP/GQ thresholds |
| PS3 | Strong | Functional studies | Not automatable — NotEvaluated |
| PS4 | Strong | Prevalence in affected vs controls | **NotEvaluated by default** — requires case-control statistics. Optional legacy proxy via `use_clinvar_stars_as_ps4_proxy` |
| PM1 | Moderate | Mutational hotspot / functional domain | ClinVar protein-position density. Capped against PP3 per Pejaver 2022 |
| PM2 | Supporting* | Absent / extremely rare in population | gnomAD raw AF, **inheritance-aware** (AD/unknown: AC=0; AR: AF ≤ 0.00007) — SVI v1.0 |
| PM3 | Supporting/Moderate/Strong/VeryStrong | In trans with pathogenic (recessive) | **SVI PM3 v1.0 points-based**: P / LP companion × phasing × hom-occurrence → 0.5 / 1.0 / 2.0 / 4.0 thresholds |
| PM4 | Moderate | Protein length change | Consequence (in-frame indel, stop-loss) |
| PM5 | Moderate | Novel missense at known pathogenic position | ClinVar protein-position index (different alt AA) |
| PM6 | Moderate | Assumed de novo (partial trio) | VCF genotype (proband + ≥1 parent). Mutually exclusive with PS2 |
| PP2 | Supporting | Missense in constrained gene | gnomAD missense Z-score ≥ 3.09 |
| PP3 | Supporting / Moderate / Strong (Pejaver 2022 + Walker 2023) | Computational pathogenic evidence | REVEL (missense only) or SpliceAI ≥ 0.2 (Supporting only) |
| PP4 | Supporting | Phenotype-specific | Not automatable — NotEvaluated |
| PP5 | Supporting | Reputable source | **Disabled by default** per ClinGen SVI |

*PM2 downgraded from Moderate to Supporting per ClinGen SVI v1.0.

#### Benign Criteria (7 automated)

| Criterion | Strength | Description | Data Source / Notes |
|-----------|----------|-------------|---------------------|
| BA1 | Standalone | Common variant (AF > 5%) | gnomAD max population AF, with **AN ≥ 2000** minimum (gnomAD v4 / SVI March 2024). Honors the **9-variant Ghosh 2018 BA1 exception list** |
| BS1 | Strong | Greater than expected frequency | gnomAD AF (gene-specific or default 0.01); same AN minimum as BA1 |
| BS2 | Strong | Observed in healthy adults | gnomAD homozygote count + OMIM inheritance |
| BS3 | Strong | Functional studies — no damage | Not automatable — NotEvaluated |
| BS4 | Strong | Lack of segregation | Not automatable — NotEvaluated |
| BP1 | Supporting | Missense in truncation-disease gene | gnomAD pLI ≥ 0.9 + misZ < 2.0 |
| BP2 | Supporting | In trans / in cis with pathogenic | Companion-variant phasing + OMIM inheritance |
| BP3 | Supporting | In-frame indel in repeat region | Consequence + RepeatMasker |
| BP4 | Supporting / Moderate / Strong / **VeryStrong** | Computational benign evidence | REVEL (missense only, **incl. VeryStrong band at REVEL ≤ 0.003**) or SpliceAI ≤ 0.1 (Walker 2023) |
| BP5 | Supporting | Alternate molecular basis | Not automatable — NotEvaluated |
| BP6 | Supporting | Reputable source — benign | **Disabled by default** per ClinGen SVI |
| BP7 | Supporting | Synonymous (mid-exon) or deep-intronic, no splice, not conserved | Consequence + SpliceAI + PhyloP + transcript exon coords. **Walker 2023**: exon-edge exclusion + deep-intronic extension |

**10 criteria require manual curation** and are marked NotEvaluated:
PS3 / PS4 (default) / BS3 / BS4 / PP1 / PP4 / PP5 (disabled) / BP2 (when
unphased) / BP5 / BP6 (disabled).

### Pejaver 2022 Calibrated REVEL Thresholds

REVEL is applied **only to missense variants** per Pejaver 2022. The
single-tool calibration replaces the previous SIFT/PolyPhen/PhyloP/GERP
ensemble (Pejaver explicitly recommends a single calibrated tool over
ad-hoc consensus).

| Direction | Strength | REVEL threshold |
|-----------|----------|-----------------|
| PP3 | Supporting | ≥ 0.644 |
| PP3 | Moderate   | ≥ 0.773 |
| PP3 | Strong     | ≥ 0.932 |
| BP4 | Supporting | ≤ 0.290 |
| BP4 | Moderate   | ≤ 0.183 |
| BP4 | Strong     | ≤ 0.016 |
| BP4 | **Very Strong** (REVEL only) | ≤ 0.003 |

A single BP4_VeryStrong is mapped to 2× `benign_strong` in the counts so
it satisfies the existing ≥2 BS → Benign rule alone (Tavtigian Bayesian
framework).

### Walker 2023 Splicing Recommendations

- **PP3 splice**: SpliceAI max_ds ≥ 0.2 → PP3 *Supporting* (no Strong from
  SpliceAI alone — Strong splice claims need experimental RNA assay).
- **BP4 splice**: SpliceAI max_ds ≤ 0.1 → BP4 Supporting.
- **Uninformative zone**: 0.10 < max_ds < 0.20 — neither fires.
- **BP7 exon-edge exclusion**: BP7 cannot fire for synonymous at first
  base or last 3 bases of an exon (canonical splice region).
- **BP7 deep-intronic extension**: BP7 may fire for intronic variants
  with offset ≥ 7 (donor side) or ≤ -21 (acceptor side) when SpliceAI is
  low and PhyloP is low.

### PVS1 Decision Tree (Abou Tayoun 2018)

| Strength | Trigger |
|----------|---------|
| **PVS1** (Very Strong) | Nonsense/frameshift predicted to undergo NMD; canonical ±1/2 splice predicted to cause NMD; whole-gene deletion in haploinsufficient gene |
| **PVS1_Strong** | NMD-escape in critical functional region |
| **PVS1_Moderate** | NMD-escape, non-critical region, ≥10% protein removed; canonical splice in last exon (NMD unlikely); start-loss with downstream Met ≤100 codons + pathogenic upstream |
| **PVS1_Supporting** | <10% protein removed in non-critical region; start-loss without strong corroborating evidence |

When NMD prediction or other transcript-level signals are missing,
PVS1 falls back to legacy full-strength VeryStrong for backward
compatibility.

### PM2 Inheritance-Aware Threshold (SVI v1.0)

| Inheritance | Threshold |
|-------------|-----------|
| AD / unknown | Strict absence (AC = 0 AND AF = 0) |
| AR | AF ≤ 0.00007 (0.007%) |
| Per-gene override | Wins over inheritance default |

Uses **raw** gnomAD AF (not FAF / popmax). Inheritance is inferred from
OMIM phenotypes.

### PM3 v1.0 Points Scoring (SVI)

| Observation | Points |
|-------------|--------|
| Confirmed in-trans + co-occurring **Pathogenic** | 1.0 |
| Confirmed in-trans + co-occurring **Likely Pathogenic** | 0.5 |
| Phase unknown + Pathogenic | 0.5 |
| Phase unknown + Likely Pathogenic | 0.25 |
| Homozygous proband observation | 0.5 (capped at 1.0 total) |

| Total points | Strength |
|--------------|----------|
| ≥ 0.5 | PM3_Supporting |
| ≥ 1.0 | PM3 (Moderate) |
| ≥ 2.0 | PM3_Strong |
| ≥ 4.0 | PM3_VeryStrong |

In-cis companions are excluded from PM3 (those count toward BP2).

### BA1 Exception List (Ghosh 2018)

Nine variants exempt from BA1 despite exceeding the 5% threshold (HFE
c.845G>A p.Cys282Tyr, GJB2 c.109G>A, F2/F5 founder alleles, etc.). Match
on `(gene_symbol, hgvs_c)`, case-insensitive. Configurable via TOML so
VCEPs can extend.

| Gene | Variant | Note |
|------|---------|------|
| ACAD9 | c.-44_-41dupTAAG | VUS |
| ACADS | c.511C>T | VUS |
| BTD | c.1330G>C | Pathogenic — biotinidase deficiency |
| GJB2 | c.109G>A | Pathogenic — DFNB1 hearing loss |
| HFE | c.187C>G | Pathogenic — hemochromatosis |
| HFE | c.845G>A | Pathogenic — hemochromatosis |
| MEFV | c.1105C>T | VUS |
| MEFV | c.1223G>A | VUS |
| PIBF1 | c.1214G>A | VUS |

### Anti-Double-Counting (PP3 Reconciliation)

A post-evaluation reconciliation pass suppresses PP3 (or PM1) under
overlap conditions called out in Pejaver 2022 and Walker 2023:

| Trigger | Suppressed | Source |
|---------|------------|--------|
| PVS1 fires AND PP3 was driven by SpliceAI | PP3 | Walker 2023 |
| PS1 fires AND PP3 was driven by REVEL | PP3 | Pejaver 2022 |
| PM5 fires AND PP3 was driven by REVEL | PP3 | Pejaver 2022 |
| PP3_Strong + PM1 (combined > Strong cap) | PM1 | Pejaver 2022 |

Suppressed criteria stay in the result with `met=false` and a
`details.suppressed_by_reconcile` note.

### gnomAD v4 AN Minimum (SVI March 2024)

BA1 and BS1 require gnomAD `all_an ≥ 2000` before they can fire. Below
the threshold the criterion is NotEvaluated rather than fired on noisy
estimates. Configurable via `min_an_for_frequency_criteria`.

### Combination Rules (19 = 18 Richards 2015 + 1 SVI Sept 2020)

**Benign:**
1. BA1 standalone → Benign
2. ≥2 BS → Benign

**Pathogenic (8):**
3. PVS + ≥1 PS
4. PVS + ≥2 PM
5. PVS + 1 PM + 1 PP
6. PVS + ≥2 PP
7. ≥2 PS
8. 1 PS + ≥3 PM
9. 1 PS + 2 PM + ≥2 PP
10. 1 PS + 1 PM + ≥4 PP

**Likely Pathogenic (7, includes ClinGen SVI Sept 2020 rule):**
11. PVS + 1 PM
12. **PVS + ≥1 PP** *(ClinGen SVI Sept 2020 — compensates PM2 downgrade; Bayesian Post_P = 0.988)*
13. 1 PS + 1–2 PM
14. 1 PS + ≥2 PP
15. ≥3 PM
16. 2 PM + ≥2 PP
17. 1 PM + ≥4 PP

**Likely Benign (2):**
18. 1 BS + 1 BP
19. ≥2 BP

**Conflict gating (PR9 fix)**: pathogenic and benign rules apply
**independently**. The result is VUS-Conflicting only when **both**
directions reach a definite call (P/LP and B/LB). Otherwise the
directional call wins.

### Trio Analysis

When a multi-sample VCF with trio configuration is provided:
- **PS2** (de novo): proband carries variant, both parents hom-ref, all
  pass DP ≥ 10 / GQ ≥ 20.
- **PM6** (assumed de novo): partial parental data; mutually exclusive
  with PS2.
- **PM3** (compound het): SVI v1.0 points scoring (above). Recessive gene
  required (OMIM).
- **BP2** (in cis/trans): for dominant genes — variant in trans with
  pathogenic; for any gene — variant in cis with pathogenic.

## ClinVar Concordance Benchmark

### Methodology

The benchmark runs `fastvep annotate --acmg --pick` end-to-end on every
ClinVar 2-star+ GRCh38 SNV / small indel and compares the issued ACMG
classification against the ClinVar review-panel call.

Concrete pipeline (`data/benchmark/run_full_benchmark.sh`):

1. **Input**: ClinVar VCF filtered to review_status ≥ 2 stars
   (`criteria_provided,_multiple_submitters,_no_conflicts` and stricter)
   on GRCh38, plus a parallel truth TSV (`chrom`, `pos`, `ref`, `alt`,
   `gene`, `clnsig`, `normalized_class`, `review_stars`, `rcv`).
2. **Annotation**: GFF3 + FASTA cache + supplementary annotation
   directory (`--sa-dir`) loaded once; all 673,660 variants annotated
   with `--acmg --pick` to a single JSON file.
3. **Concordance** (`analysis/acmg_benchmark/real_data/03_evaluate_concordance.py`):
   stream-parses the JSON via `ijson` (memory-bounded — output is ~24 GB
   pretty-printed), keys each variant on `(chrom, pos, ref, alt)`,
   reads the picked transcript's `acmg.classification`, and fills a
   5×6 truth × predicted matrix (extra column for NoCall).

Outputs: `concordance_summary.txt` (free-text rollup),
`concordance_matrix.csv`, `concordance_by_chrom.csv`,
`concordance_by_consequence.csv`, `criterion_firing_rates.csv`,
`rule_distribution.csv`, `discrepancies.tsv` (top 10k
opposite-direction calls).

### Supplementary Annotation Build

| Source | Build | Records |
|--------|-------|---------|
| ClinVar (.osa) | `fastvep sa-build --source clinvar` from `clinvar.vcf.gz` | 4,402,501 |
| ClinVar protein (.oga) | `--source clinvar_protein` from `variant_summary.txt.gz` | 4,554 genes |
| gnomAD v4.1 exomes (.osa, per-chrom) | tabix-extracted to ClinVar 2-star+ regions (24,350 merged ranges, `bedtools merge -d 5000`), `--source gnomad` per chrom (chr1..22, X, Y, MT) | 25 × .osa |
| gnomAD v4.1 gene constraints (.oga) | `--source gnomad_gene` from `gnomad.v4.1.constraint_metrics.tsv` | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom) | per-chromosome split from `revel-v1.3_all_chromosomes.zip` to bound RAM | 24 × .osa |

The gnomAD bulk-extraction path uses `tabix` against the public bgz on
`gs://gcp-public-data--gnomad/...`. We tested the gnomAD GraphQL API as
an alternative: it is fine for ad-hoc per-variant lookups, but
**rate-limits aggressively (HTTP 429) even single-threaded with 5-try
exponential backoff**, so it cannot replace tabix for the 24 K-region
extraction.

### Speed (single host, Apple Silicon, release build)

| Stage | Time | Throughput |
|-------|------|-----------|
| `fastvep annotate --acmg --pick` on 673,660 variants | **2,591 s (43 min)** | **260 variants/s** |
| Streaming concordance parse of 24 GB JSON via ijson | ~12 min | — |

(All 25 SA databases loaded once at process start; 99 % CPU during the
annotation phase.)

### Real-Data Concordance (ClinVar 2-star+, April 2026 release)

Truth records: **673,660** · Classified: **627,375** · Truth-not-annotated: **46,285** (variants where `--pick` selected a transcript without an ACMG block — typically intergenic / regulatory regions where no canonical-transcript context exists).

#### Truth × predicted matrix

| Truth ↓ / Predicted → | P | LP | VUS | LB | B | NoCall |
|--|--:|--:|--:|--:|--:|--:|
| Pathogenic (n=49,882) | 80 | 7,756 | 42,025 | 12 | 8 | 1 |
| Likely Pathogenic (n=11,589) | 5 | 2,418 | 9,161 | 5 | 0 | 0 |
| VUS (n=288,912) | 1 | 156 | 279,031 | 7,749 | 1,975 | 0 |
| Likely Benign (n=126,036) | 0 | 1 | 121,970 | 778 | 3,287 | 0 |
| Benign (n=150,957) | 0 | 5 | 100,833 | 1,837 | 48,282 | 0 |

#### Headline metrics

| Metric | Value |
|--------|------:|
| Exact-match (truth = predicted) | 52.7 % |
| Same-direction (truth & predicted both P-tier or both B-tier or both VUS) | 54.7 % |
| Opposite-direction (P/LP truth → B/LB predicted, or vice versa) | **31 / 627,375 = 0.005 %** |
| NoCall after annotation | 0.0 % |

Per-class same-direction recall:

| Truth class | Same-dir count | Recall |
|-------------|---------------:|------:|
| Pathogenic | 7,836 / 49,882 | **15.7 %** |
| Likely Pathogenic | 2,423 / 11,589 | **20.9 %** |
| VUS | 279,031 / 288,912 | **96.6 %** |
| Likely Benign | 4,065 / 126,036 | **3.2 %** |
| Benign | 50,119 / 150,957 | **33.2 %** |

#### Most-triggered combination rules

| Rule | Count |
|------|------:|
| BA1 | 42,081 |
| ≥2 BS | 11,471 |
| PS + 1–2 PM | 9,463 |
| BS + BP | 7,713 |
| ≥2 BP | 2,668 |
| PVS + PM | 503 |
| PVS + ≥1 PP (SVI Sept 2020) | 304 |
| PS + 2 PM + ≥2 PP | 62 |

### Interpretation

- **High specificity, low sensitivity for P/LP** is expected and by
  design. The classifier is automation-only; the criteria carrying the
  highest evidence weight in real curation — PS3 (functional), BS3
  (functional non-damaging), BS4 (lack of segregation), PP1 (segregation),
  PP4 (phenotype-specific), PP5/BP6 (reputable source, disabled per
  ClinGen SVI) — are unavailable from variant-level data alone, so most
  P/LP variants drop to VUS for lack of those signals.
- **Opposite-direction rate is 31 / 627,375 (0.005 %)**: when the
  classifier *does* commit to a directional call, it almost never
  contradicts the curated review-panel call. Discrepancies are listed
  in `data/benchmark/output_full/discrepancies.tsv` for case-by-case
  review.
- **Likely_benign collapses into VUS** (3.2 % recall): without PS3/BS3
  data, BP4 is the main driver toward benign — and when only one
  benign-supporting criterion fires it is sub-threshold for any benign
  rule (≥1 BS + ≥1 BP, or ≥2 BP), so the call falls to VUS.
- **Pathogenic exact-match (80 / 49,882) vs. same-direction (7,836 /
  49,882)** shows that of the variants the classifier reaches P/LP on,
  it lands at LP rather than P most of the time — consistent with the
  ClinGen SVI Sept 2020 PVS+PP rule which the classifier triggers
  heavily (304 firings) but which only escalates to LP, not P.

### Two diagnostic findings worth flagging

The matrix is dominated by two patterns that look like classifier bugs
but aren't — they're both downstream of missing SA data sources rather
than incorrect criterion logic. Documenting them here so future
benchmark deltas are interpretable.

**(1) Likely-benign collapse to VUS (3.2 % LB recall)**

Most LB truth is synonymous (62 K) or intronic (25 K). For these
classes, BP7 is the canonical benign-supporting criterion. BP7 fires
**zero times** in this run because Walker 2023 requires SpliceAI ≤ 0.2
*and* PhyloP < 2.0, and **neither SpliceAI nor PhyloP SA databases were
loaded**. With BP7 silenced, synonymous LB has no path to the ≥2 BP
or BS+BP rules; everything falls to VUS. (BP4 missense / REVEL fires
on the missense subset of LB only.)

**(2) Pathogenic exact-match (80 / 49,882) vs same-direction (7,836 / 49,882)**

PVS1 fires only on 10.5 % of pathogenic-truth variants because
`is_lof_intolerant_gene` requires gnomAD `pLI ≥ 0.9` *or* `LOEUF ≤
0.35` *or* an OMIM phenotype, and **OMIM was not loaded** in this run.
Many stop_gained / frameshift pathogenic ClinVar entries are in genes
that don't meet the strict gnomAD thresholds (rare LOF tolerance,
small-N constraint), so they get filtered out without OMIM as a
backstop. Where PVS1 *does* fire, the partner criterion gating is
PM2_Supporting (9.0 % on Pathogenic) — the rate-limiting step for the
PVS+PP rule (304 firings) which lands at LP rather than P.

### Data sources that would change the picture if loaded

| SA source | Affects | Expected delta |
|-----------|---------|----------------|
| **PhyloP / GERP** (.osa) | BP7 firing on synonymous/deep-intronic | Should lift LB recall from 3.2 % toward ~30 % |
| **SpliceAI** (.osa) | PP3 splice, BP4 splice, BP7 splice gate, PVS1 splice grading | Lifts P recall on intronic/canonical-splice; flips many intronic VUS → P/LP or B |
| **OMIM** (.oga) | PVS1 LOF-gene fallback (alongside pLI/LOEUF) | Lifts P recall on stop_gained / frameshift in non-constrained disease genes |

These are infrastructure (data) rather than code changes — the
classifier already consumes all four when available (`phylop`, `gerp`,
`spliceai`, `omim` json keys in `sa_extract.rs`). We simply did not
build them for this benchmark.

### Limitations of the automated benchmark

1. **Inherently conservative**: PS3/BS3/BS4/PP1/PP4/BP5 are all
   NotEvaluated. Manual curators outperform any variant-level automation
   for these categories by design. The benchmark measures
   classifier-vs-curator agreement, not classifier-vs-ground-truth.
2. **PVS1 / PS1 / BP7 fallback paths**: when the pipeline cannot
   compute Abou Tayoun decision-tree signals (NMD, %protein removed)
   for a specific transcript, those criteria fall back to conservative
   legacy behavior. PVS1_Strong / PVS1_Moderate / PVS1_Supporting
   firings in the table reflect cases where the pipeline *did* derive
   the tree signal.
3. **PS4 NotEvaluated by default**: the previous ClinVar-stars proxy was
   replaced; opt back in via `use_clinvar_stars_as_ps4_proxy` for a
   backward-comparable benchmark.
4. **gnomAD v4 mid / remaining populations**: added to the parser and
   `max_pop_af` after the audit. The 9 chromosome `.osa` files built
   before this change (chr 6, 13, 18, 20, 21, 22, MT, X, Y) lack those
   keys; impact is small (mid + remaining ≈ 5 % of v4 cohort).

## Configuration

```toml
# Frequency thresholds
ba1_af_threshold = 0.05
bs1_af_threshold = 0.01
pm2_af_threshold = 0.0001            # legacy single-threshold field (back-compat)
pm2_ad_af_threshold = 0.0            # AD / unknown: strict absence
pm2_ar_af_threshold = 0.00007        # AR threshold (SVI v1.0)
min_an_for_frequency_criteria = 2000 # gnomAD v4 AN minimum (SVI March 2024)

# REVEL thresholds (Pejaver 2022; missense only)
pp3_revel_supporting = 0.644
pp3_revel_moderate = 0.773
pp3_revel_strong = 0.932
bp4_revel_supporting = 0.290
bp4_revel_moderate = 0.183
bp4_revel_strong = 0.016
bp4_revel_very_strong = 0.003        # only REVEL reaches this band

# SpliceAI thresholds (Walker 2023)
spliceai_pathogenic = 0.2
spliceai_benign = 0.1

# Conservation
phylop_conserved = 2.0

# Gene constraints
pli_lof_intolerant = 0.9
loeuf_lof_intolerant = 0.35
pp2_misz_threshold = 3.09
pm1_hotspot_window = 5
pm1_hotspot_min_pathogenic = 3

# ClinGen SVI behavior modifications
pm2_downgrade_to_supporting = true
use_pp5_bp6 = false
use_clinvar_stars_as_ps4_proxy = false

# BA1 exception list — defaults to the 9-variant Ghosh 2018 set;
# users can extend or replace via TOML.
[[ba1_exceptions]]
gene = "HFE"
hgvs_c = "c.845G>A"
reason = "Hereditary hemochromatosis"

# Gene-specific overrides
[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001

# Per-disorder overrides for multi-disorder genes (SVI July 2025 scaffold)
[gene_overrides.GENE_X.disorders.disorder_a]
inheritance = "AR"
pm2_af_threshold = 0.00007
```

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424.
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat*. 2018;39(11):1517-1524.
- Ghosh R, et al. Updated recommendation for the benign stand-alone ACMG/AMP criterion. *Hum Mutat*. 2018;39(11):1525-1530.
- ClinGen SVI Recommendation for Absence/Rarity (PM2) — Version 1.0. Approved September 4, 2020.
- ClinGen SVI Recommendation for In-Trans Criterion (PM3) — Version 1.0.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163-2177.
- Walker LC, et al. (ClinGen SVI Splicing Subgroup). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing. *Am J Hum Genet*. 2023;110(7):1046-1067.
- ClinGen SVI Working Group. Guidance to VCEPs Regarding gnomAD v4 (March 2024).
- ClinGen SVI Working Group. Guidance Classifying Variants in Genes Associated with Multiple Disorders (July 2025).
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med*. 2018;20(9):1054-1060.
- Lek M, et al. Analysis of protein-coding genetic variation in 60,706 humans. *Nature*. 2016;536(7616):285-291.

# Benchmark setup state

This file is **gitignored** — it tracks what's been downloaded / built locally so we don't re-do completed steps. Update it as you complete each step.

## Genome references (`refs/`)

Reuses files already on disk under `test_data/organisms/human/`:

| File | Path | Size |
|------|------|------|
| GFF3 | test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3 | 1.22 GB |
| GFF3 cache | test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3.fastvep.cache | 388 MB |
| FASTA | test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa | 3.15 GB |
| FASTA index | test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai | 6 KB |

## SA sources (`sa_sources/`)

| File | Status | Size | Used by |
|------|--------|------|---------|
| clinvar.vcf.gz (from `test_data/sa_databases/`) | ✅ | 190 MB | clinvar.osa |
| variant_summary.txt.gz | ✅ | 419 MB | clinvar_protein.oga |
| gnomad.v4.1.constraint_metrics.tsv | ✅ | 95 MB | gnomad_genes.oga |
| revel-v1.3_all_chromosomes.zip | ✅ | 637 MB | per-chrom revel.osa |
| revel_per_chrom/revel.chr*.csv | ✅ | 6.5 GB | per-chrom revel.osa |
| clinvar_2star.vcf | ✅ | 394 MB | benchmark input |
| clinvar_2star_truth.tsv | ✅ | 47 MB | benchmark truth |
| clinvar_2star_regions.bed | ✅ | 800 KB | gnomAD tabix extraction (24350 ranges) |
| gnomad_extracts/gnomad_chr*.vcf | 🔧 in progress | varies | per-chrom gnomad.osa |

## Built SA databases (`sa_db/`)

| Source | Status | Records | Path |
|--------|--------|---------|------|
| clinvar.osa | ✅ | 4,402,501 | sa_db/clinvar.osa |
| clinvar_protein.oga | ✅ | 4,554 | sa_db/clinvar_protein.oga |
| gnomad_genes.oga | ✅ | 18,173 | sa_db/gnomad_genes.oga |
| revel_chr1..22,X,Y.osa | ✅ | 24 chrom files | sa_db/revel_chr*.osa |
| gnomad_chr22.osa | 🔧 in progress | — | sa_db/gnomad_chr22.osa |
| gnomad_chr1..21,X,Y.osa | ⏳ pending | — | — |

## Benchmark runs

| Run | Date | SA loaded | Outcome | Output |
|-----|------|-----------|---------|--------|
| chr17 interim | 2026-04-29 | clinvar+protein+genes+revel(none yet) | partial — see analysis/acmg_benchmark/real_data/output_chr17/ | committed |
| (next) full | TBD | clinvar+protein+genes+revel(all)+gnomad(extracted) | TBD | data/benchmark/output_v7/ |

## Convention

- ✅ done
- 🔧 in progress
- ⏳ pending
- ❌ broken / needs redo

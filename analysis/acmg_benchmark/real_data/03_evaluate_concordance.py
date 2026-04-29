#!/usr/bin/env python3
"""
Evaluate ACMG classifier concordance against ClinVar 2-star+ truth.

Inputs
  --truth        TSV: chrom\tpos\tref\talt\tgene\tclnsig\tnormalized_class\treview_stars\trcv
  --predictions  fastvep --output-format json (single JSON array; each
                 element has transcript_consequences[].acmg)

Outputs (under --out)
  concordance_matrix.csv             5×6 (truth × predicted+NoCall)
  concordance_by_chrom.csv           per-chromosome concordance
  concordance_by_consequence.csv     top consequences × class concordance
  concordance_summary.txt            full text report
  criterion_firing_rates.csv         per-criterion fire rates by truth class
  rule_distribution.csv              triggered-rule frequencies
  discrepancies.tsv                  opposite-direction calls (top 10k)
"""

import argparse
import json
import csv
from pathlib import Path
from collections import defaultdict, Counter

import ijson

CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]


def load_truth(path):
    truth = {}
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            truth[key] = row
    return truth


def variant_keys(rec):
    chrom = str(rec.get("seq_region_name") or "")
    pos = str(rec.get("start") or rec.get("position") or "")
    allele = rec.get("allele_string", "") or ""
    if "/" not in allele:
        return []
    ref, alts = allele.split("/", 1)
    return [(chrom, pos, ref, a) for a in alts.split(",")]


def class_label(c):
    if c == "Pathogenic":
        return "Pathogenic"
    if c in ("Likely_pathogenic", "LikelyPathogenic"):
        return "Likely_pathogenic"
    if c in ("Uncertain_significance", "UncertainSignificance", "VUS"):
        return "VUS"
    if c in ("Likely_benign", "LikelyBenign"):
        return "Likely_benign"
    if c == "Benign":
        return "Benign"
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--truth", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--out", default="output")
    args = p.parse_args()

    truth = load_truth(args.truth)
    print(f"Loaded {len(truth):,} truth records")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    cm = {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES}
    cm_consequence = defaultdict(lambda: {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES})
    cm_chrom = defaultdict(lambda: {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES})
    criterion_fires = {tc: Counter() for tc in CLASSES}
    criterion_evaluable = {tc: Counter() for tc in CLASSES}
    rule_dist = Counter()
    discrepancies = []
    matched = set()
    n_classified = 0

    # Stream-parse: fastvep --output-format json emits one JSON array
    # (pretty-printed) which can run into 20+ GB on a full ClinVar 2-star+
    # set. ijson keeps memory bounded; line-delimited JSONL also works.
    def _records(path):
        with open(path) as fh:
            first = fh.read(1)
            fh.seek(0)
            if first == "[":
                yield from ijson.items(fh, "item")
            else:
                for line in fh:
                    line = line.strip()
                    if line:
                        yield json.loads(line)

    for rec in _records(args.predictions):
        for k in variant_keys(rec):
            if k not in truth:
                continue
            t = truth[k]
            tc_truth = t["normalized_class"]
            acmg = None
            top_csq = "unknown"
            for tc in rec.get("transcript_consequences", []) or []:
                if "acmg" in tc:
                    acmg = tc["acmg"]
                    cs = tc.get("consequence_terms", [])
                    if cs:
                        top_csq = cs[0]
                    break
            if acmg is None:
                cm[tc_truth]["NoCall"] += 1
                cm_chrom[k[0]][tc_truth]["NoCall"] += 1
                continue
            pc = class_label(acmg.get("classification")) or "NoCall"
            cm[tc_truth][pc] += 1
            cm_chrom[k[0]][tc_truth][pc] += 1
            cm_consequence[top_csq][tc_truth][pc] += 1
            rule = acmg.get("triggered_rule") or ""
            if rule:
                rule_dist[rule] += 1
            for c in acmg.get("criteria", []) or []:
                if c.get("evaluated"):
                    criterion_evaluable[tc_truth][c["code"]] += 1
                if c.get("met"):
                    criterion_fires[tc_truth][c["code"]] += 1
            n_classified += 1
            matched.add(k)
            if (
                tc_truth in ("Pathogenic", "Likely_pathogenic")
                and pc in ("Benign", "Likely_benign", "VUS", "NoCall")
            ) or (
                tc_truth in ("Benign", "Likely_benign")
                and pc in ("Pathogenic", "Likely_pathogenic")
            ):
                if len(discrepancies) < 10000:
                    met_codes = ";".join(c["code"] for c in (acmg.get("criteria") or []) if c.get("met"))
                    discrepancies.append((
                        k[0], k[1], k[2], k[3], t["gene"], t["review_stars"],
                        tc_truth, pc, top_csq, rule, met_codes
                    ))
            break

    n_truth = len(truth)
    n_unmatched = n_truth - len(matched)

    # ── concordance_matrix.csv ──
    matrix_path = out_dir / "concordance_matrix.csv"
    with matrix_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["truth"] + CLASSES + ["NoCall"])
        for tcl in CLASSES:
            w.writerow([tcl] + [cm[tcl][pc] for pc in CLASSES + ["NoCall"]])

    # ── concordance_by_chrom.csv ──
    by_chrom_path = out_dir / "concordance_by_chrom.csv"
    with by_chrom_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["chrom", "truth", "n", "exact", "same_dir", "opposite", "no_call"])
        for chrom in sorted(cm_chrom.keys(), key=lambda x: (len(x), x)):
            for tcl in CLASSES:
                row = cm_chrom[chrom][tcl]
                n = sum(row.values())
                if n == 0:
                    continue
                exact = row[tcl]
                if tcl in ("Pathogenic", "Likely_pathogenic"):
                    same = row["Pathogenic"] + row["Likely_pathogenic"]
                    opp = row["Benign"] + row["Likely_benign"]
                elif tcl == "VUS":
                    same = row["VUS"]
                    opp = 0
                else:
                    same = row["Benign"] + row["Likely_benign"]
                    opp = row["Pathogenic"] + row["Likely_pathogenic"]
                w.writerow([chrom, tcl, n, exact, same, opp, row["NoCall"]])

    # ── concordance_by_consequence.csv ──
    consq_counts = sorted(
        cm_consequence.items(),
        key=lambda kv: -sum(sum(d.values()) for d in kv[1].values()),
    )[:15]
    by_csq_path = out_dir / "concordance_by_consequence.csv"
    with by_csq_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["consequence", "truth"] + CLASSES + ["NoCall", "n"])
        for csq, mat in consq_counts:
            for tcl in CLASSES:
                row = mat[tcl]
                n = sum(row.values())
                if n == 0:
                    continue
                w.writerow([csq, tcl] + [row[pc] for pc in CLASSES + ["NoCall"]] + [n])

    # ── criterion_firing_rates.csv ──
    fr_path = out_dir / "criterion_firing_rates.csv"
    all_codes = sorted({c for cnt in criterion_fires.values() for c in cnt} |
                       {c for cnt in criterion_evaluable.values() for c in cnt})
    with fr_path.open("w") as f:
        w = csv.writer(f)
        header = ["criterion"]
        for tcl in CLASSES:
            header.extend([f"{tcl}_evaluable", f"{tcl}_fired", f"{tcl}_pct"])
        w.writerow(header)
        for code in all_codes:
            r = [code]
            for tcl in CLASSES:
                ev = criterion_evaluable[tcl].get(code, 0)
                fi = criterion_fires[tcl].get(code, 0)
                pct = (fi / ev * 100) if ev else 0
                r.extend([ev, fi, f"{pct:.1f}"])
            w.writerow(r)

    # ── rule_distribution.csv ──
    rd_path = out_dir / "rule_distribution.csv"
    with rd_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["rule", "n"])
        for rule, n in rule_dist.most_common():
            w.writerow([rule, n])

    # ── discrepancies.tsv ──
    disc_path = out_dir / "discrepancies.tsv"
    with disc_path.open("w") as f:
        f.write("chrom\tpos\tref\talt\tgene\tstars\ttruth\tpredicted\tconsequence\trule\tmet_criteria\n")
        for d in discrepancies:
            f.write("\t".join(str(x) for x in d) + "\n")

    # ── concordance_summary.txt ──
    totals = {"n": 0, "exact": 0, "same": 0, "opp": 0, "nc": 0}
    summary_path = out_dir / "concordance_summary.txt"
    with summary_path.open("w") as f:
        f.write("ClinVar 2-star+ concordance against fastvep ACMG classifier (real data)\n")
        f.write("=" * 75 + "\n\n")
        f.write(f"Truth records:       {n_truth:,}\n")
        f.write(f"Classified:          {n_classified:,}\n")
        f.write(f"Truth not annotated: {n_unmatched:,}\n\n")
        f.write("Per-class breakdown (entire dataset):\n")
        f.write(f"{'truth':<22} {'n':>8} {'exact':>8} {'same_dir':>10} {'opposite':>10} {'no_call':>8}\n")
        for tcl in CLASSES:
            row = cm[tcl]
            n = sum(row.values())
            exact = row[tcl]
            if tcl in ("Pathogenic", "Likely_pathogenic"):
                same = row["Pathogenic"] + row["Likely_pathogenic"]
                opp = row["Benign"] + row["Likely_benign"]
            elif tcl == "VUS":
                same = row["VUS"]
                opp = 0
            else:
                same = row["Benign"] + row["Likely_benign"]
                opp = row["Pathogenic"] + row["Likely_pathogenic"]
            no_call = row["NoCall"]
            f.write(f"{tcl:<22} {n:>8} {exact:>8} {same:>10} {opp:>10} {no_call:>8}\n")
            for k, v in [("n", n), ("exact", exact), ("same", same), ("opp", opp), ("nc", no_call)]:
                totals[k] += v
        f.write(
            f"\n{'TOTAL':<22} {totals['n']:>8} {totals['exact']:>8} "
            f"{totals['same']:>10} {totals['opp']:>10} {totals['nc']:>8}\n"
        )
        if totals["n"]:
            f.write(f"\nExact-match rate:        {totals['exact']/totals['n']*100:.1f}%\n")
            f.write(f"Same-direction rate:     {totals['same']/totals['n']*100:.1f}%\n")
            f.write(f"Opposite-direction rate: {totals['opp']/totals['n']*100:.1f}%\n")
            f.write(f"NoCall rate:             {totals['nc']/totals['n']*100:.1f}%\n")

        f.write("\nTop 25 triggered rules:\n")
        for rule, n in rule_dist.most_common(25):
            f.write(f"  {n:>8}  {rule}\n")

        f.write("\nCriterion firing rates (% of variants where the criterion was evaluable AND fired):\n")
        f.write(f"{'criterion':<22}")
        for tcl in CLASSES:
            f.write(f"  {tcl:>22}")
        f.write("\n")
        for code in all_codes:
            f.write(f"{code:<22}")
            for tcl in CLASSES:
                ev = criterion_evaluable[tcl].get(code, 0)
                fi = criterion_fires[tcl].get(code, 0)
                pct = (fi / ev * 100) if ev else 0
                f.write(f"  {fi:>6}/{ev:>6} ({pct:>5.1f}%)")
            f.write("\n")

    print("\nOutputs:")
    for p in (matrix_path, by_chrom_path, by_csq_path, fr_path, rd_path, disc_path, summary_path):
        print(f"  {p}")
    print()
    print(open(summary_path).read())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create manuscript rate-table artifacts from the exported CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

HEADERS = [
    "Arm",
    "Metric",
    "Unit",
    "n wells",
    "n elec.",
    "Baseline median",
    "Treatment median",
    "Baseline mean +/- SEM",
    "Treatment mean +/- SEM",
    "Median change (%)",
    "p (well)",
]


CAPTION = (
    "Summary of firing and burst metrics. Descriptive values "
    "are pooled electrode-level summaries using the same inclusion rules as "
    "Figures 2 and 4. DOI p values are exact paired Wilcoxon signed-rank "
    "tests on per-well median percent change. The ketanserin arm is "
    "reported descriptively because n = 3 wells, so no p values are shown."
)


def fmt(value: str) -> str:
    x = float(value)
    if abs(x) >= 100:
        return f"{x:.1f}"
    if abs(x) >= 10:
        return f"{x:.2f}"
    return f"{x:.3f}"


def fmt_p(value: str) -> str:
    x = float(value)
    if x < 0.001:
        return "<0.001"
    return f"{x:.3f}"


def read_rows(csv_path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with csv_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            p_value = (
                "" if r["study"] == "ket" else fmt_p(r["p_wilcoxon_well_pct_change"])
            )
            rows.append(
                [
                    r["treatment"],
                    r["metric"],
                    r["unit"],
                    r["n_wells"],
                    r["n_electrode_pairs"],
                    fmt(r["baseline_median"]),
                    fmt(r["treatment_median"]),
                    f'{fmt(r["baseline_mean"])} +/- {fmt(r["baseline_sem"])}',
                    f'{fmt(r["treatment_mean"])} +/- {fmt(r["treatment_sem"])}',
                    fmt(r["per_well_median_pct_change"]),
                    p_value,
                ]
            )
    return rows


def latex_escape(text: str) -> str:
    return (
        text.replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("_", r"\_")
        .replace("#", r"\#")
        .replace("{", r"\{")
        .replace("}", r"\}")
        .replace("<", r"$<$")
        .replace(">", r"$>$")
        .replace("+/-", r"$\pm$")
    )


def write_latex(rows: list[list[str]], tex_path: Path) -> None:
    tex_path.parent.mkdir(parents=True, exist_ok=True)
    colspec = "llrrrrrrrrr"
    with tex_path.open("w", encoding="utf-8") as fh:
        fh.write(
            "\\documentclass[10pt]{article}\n"
            "\\usepackage[margin=0.45in,landscape]{geometry}\n"
            "\\usepackage{booktabs}\n"
            "\\usepackage{array}\n"
            "\\usepackage{caption}\n"
            "\\captionsetup{font=small,labelfont=bf}\n"
            "\\begin{document}\n"
            "\\pagestyle{empty}\n"
            "\\begin{table}[ht]\n"
            "\\centering\n"
            "\\caption{Summary of firing and burst metrics.}\n"
            "\\scriptsize\n"
            f"\\begin{{tabular}}{{{colspec}}}\n"
            "\\toprule\n"
        )
        fh.write(" & ".join(latex_escape(h) for h in HEADERS) + " \\\\\n")
        fh.write("\\midrule\n")
        for row in rows:
            fh.write(" & ".join(latex_escape(v) for v in row) + " \\\\\n")
        fh.write(
            "\\bottomrule\n"
            "\\end{tabular}\n"
            "\\caption*{Descriptive values are pooled electrode-level summaries using the same inclusion rules as Figures 2 and 4. DOI p values are exact paired Wilcoxon signed-rank tests on per-well median percent change. The ketanserin arm is reported descriptively because n = 3 wells, so no p values are shown.}\n"
            "\\end{table}\n"
            "\\end{document}\n"
        )


def update_spec(spec_path: Path, rows: list[list[str]]) -> None:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    content = spec["content"]
    content[:] = [item for item in content if item.get("id") != "table_rate_results"]
    table = {
        "type": "table",
        "id": "table_rate_results",
        "caption": CAPTION,
        "headers": HEADERS,
        "rows": rows,
    }
    insert_after = next(
        (i for i, item in enumerate(content) if item.get("id") == "fig4"),
        len(content) - 1,
    )
    content.insert(insert_after + 1, table)
    spec_path.write_text(
        json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="output/tables/rate_results_table.csv")
    parser.add_argument("--tex", default="paper/tables/rate_results_table.tex")
    parser.add_argument("--spec", default="paper/spec.json")
    parser.add_argument("--update-spec", action="store_true")
    args = parser.parse_args()

    rows = read_rows(Path(args.csv))
    write_latex(rows, Path(args.tex))
    if args.update_spec:
        update_spec(Path(args.spec), rows)


if __name__ == "__main__":
    main()

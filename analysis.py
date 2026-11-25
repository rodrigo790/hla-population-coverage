"""Single-locus carrier coverage in published NMDP donor cohorts."""
import numpy as np
import pandas as pd

from common import RESULTS, download, plot_style, save_summary

PUBLICATION = (
    "https://network.nmdp.org/services-support/bioinformatics-immunobiology/"
    "haplotype-frequencies/high-resolution-hla-alleles-and-haplotypes-in-the-us-population"
)
SOURCES = {"hla.xls": (
    "https://network.nmdp.org/-/media/project/nmdp/networkresourcecenter/"
    "bioinformatics/a-xls.xls?rev=df43b341acee4d47908ed105e19a52d8"
    "&hash=898B1B5EC5A7376854E29A02BA7E7065"
)}
COHORTS = ("EUR", "AFA", "API", "HIS")
PANELS = {
    "A*02:01g": ("0201g",),
    "A*01:01g + A*02:01g": ("0101g", "0201g"),
    "A*01:01g + A*02:01g + A*24:02g": ("0101g", "0201g", "2402g"),
}


def carrier_fraction(p, f=0):
    """Probability of at least one selected haplotype in a diploid single locus."""
    if not 0 <= p <= 1 or not 0 <= f <= 1:
        raise ValueError("p and f must be fractions between zero and one")
    return 1 - (1-p)**2 - f*p*(1-p)


def read_frequencies(path):
    table = pd.read_excel(path, sheet_name="A", dtype={"A A": str})
    allele_column = table.columns[0]
    table[allele_column] = table[allele_column].astype(str).str.strip()
    frequencies = table.set_index(allele_column)[[f"{group}_freq" for group in COHORTS]]
    frequencies.columns = list(COHORTS)
    frequencies = frequencies.apply(pd.to_numeric, errors="coerce")
    return frequencies


def scenarios(frequencies):
    rows = []
    for name, panel in PANELS.items():
        for group in COHORTS:
            selected = frequencies.loc[list(panel), group]
            if len(selected) != len(panel) or selected.isna().any():
                raise ValueError(f"Missing allele group from {group}")
            mass = float(selected.sum())
            rows.extend({
                "panel": name, "donor_group": group, "F": f, "allele_mass": mass,
                "carrier_coverage": carrier_fraction(mass, f)
            } for f in (0, 0.05, 0.10))
    return pd.DataFrame(rows)


def main():
    RESULTS.mkdir(exist_ok=True)
    path = download("hla.xls", SOURCES["hla.xls"])
    table = scenarios(read_frequencies(path))
    table.to_csv(RESULTS / "scenarios.csv", index=False, float_format="%.6f")
    base = table[table.F.eq(0)].pivot(index="panel", columns="donor_group", values="carrier_coverage")
    base = base.reindex(index=PANELS.keys(), columns=COHORTS)
    plt = plot_style()
    fig, ax = plt.subplots(figsize=(8.7, 3.4), constrained_layout=True)
    heat = ax.imshow(base, cmap="GnBu", vmin=0, vmax=1)
    for i, row in enumerate(base.to_numpy()):
        for j, value in enumerate(row):
            ax.text(j, i, f"{value:.1%}", ha="center", va="center",
                    color="#f7fafc" if value > 0.65 else "#173049")
    ax.set(xticks=range(len(COHORTS)), xticklabels=COHORTS, yticks=range(len(PANELS)),
           yticklabels=PANELS.keys(),
           title="Expected HLA-A carrier coverage · Hardy–Weinberg")
    ax.grid(False)
    fig.colorbar(heat, ax=ax, label="Carrier fraction")
    fig.savefig(RESULTS / "coverage.png")
    plt.close(fig)
    widest = base.max(axis=1) - base.min(axis=1)
    save_summary({
        "cohorts": len(COHORTS), "panels": len(PANELS), "scenarios": len(table),
        "largest_between_group_gap_pp": round(float(widest.max() * 100), 2),
        "coverage_three_group_panel": {group: round(float(value), 4)
                                       for group, value in base.iloc[-1].items()},
    })


if __name__ == "__main__":
    main()

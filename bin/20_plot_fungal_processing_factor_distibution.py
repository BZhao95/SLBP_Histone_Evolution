#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
from matplotlib.patches import Rectangle


DEFAULT_CLADES = [
    "Basidiomycota",
    "Ascomycota",
    "Mucoromycota",
    "Zoopagomycota",
    "Blastocladiomycota",
    "Chytridiomycota",
    "Cryptomycota",
    "Microsporidia",
]

FACTOR_ORDER = [
    "SLBP",
    "LSM10",
    "LSM11",
    "U7-like",
    "GLD2-like",
    "Trf4-like",
]

TEXT_TAXID_FILES = {
    "SLBP": Path("slbp/slbp_taxids.txt"),
    "LSM10": Path("lsm10/lsm10_taxids.txt"),
    "LSM11": Path("lsm11/lsm11_taxids.txt"),
    "U7-like": Path("GLD2/u7_taxids.txt"),
}

TSV_TAXID_FILES = {
    "GLD2-like": Path("GLD2/gld-2_A-adding_accessions_species_protein_clades.tsv"),
    "Trf4-like": Path("GLD2/trf4_accessions_species_protein_clades.tsv"),
}

COLOR_MAP = mpl.colors.LinearSegmentedColormap.from_list(
    "fungal_factor_abundance",
    ["#fff4df", "#ffd591", "#ffb347", "#ed8b00", "#b85c00", "#783400"],
)


def get_args():
    parser = argparse.ArgumentParser(
        description="Count fungal histone-processing factors and plot their distribution."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing the slbp, lsm10, lsm11 and GLD2 subdirectories.",
    )
    parser.add_argument(
        "-c",
        "--counts-output",
        type=Path,
        default=Path("fungal_processing_factor_counts.csv"),
        help="Output CSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("fungal_processing_factor_distribution.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--clades",
        nargs="+",
        default=DEFAULT_CLADES,
        help="Clades in top-to-bottom plotting order.",
    )
    parser.add_argument(
        "--color-max",
        type=float,
        default=100,
        help="Count mapped to the darkest colour.",
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--fig-width", type=float, default=13)
    parser.add_argument("--fig-height", type=float, default=7.5)
    return parser.parse_args()


def load_taxids_from_text(file_path):
    if not file_path.is_file():
        raise FileNotFoundError(f"Taxid file not found: {file_path}")

    taxids = set()

    with file_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value:
                continue

            try:
                taxids.add(int(value))
            except ValueError:
                print(
                    f"Warning: ignoring invalid taxid '{value}' in "
                    f"{file_path}, line {line_number}"
                )

    return taxids


def load_taxids_from_third_column(file_path):
    if not file_path.is_file():
        raise FileNotFoundError(f"TSV file not found: {file_path}")

    table = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        dtype=str,
        comment="#",
    )

    if table.shape[1] < 3:
        raise ValueError(
            f"{file_path} has {table.shape[1]} columns; at least three are required."
        )

    taxids = pd.to_numeric(table.iloc[:, 2], errors="coerce").dropna()
    return set(taxids.astype(int))


def load_factor_taxids(input_dir):
    factor_taxids = {
        factor: load_taxids_from_text(input_dir / relative_path)
        for factor, relative_path in TEXT_TAXID_FILES.items()
    }
    factor_taxids.update(
        {
            factor: load_taxids_from_third_column(input_dir / relative_path)
            for factor, relative_path in TSV_TAXID_FILES.items()
        }
    )
    return factor_taxids


def get_clade_taxids(ncbi, clade):
    translation = ncbi.get_name_translator([clade])

    if clade not in translation:
        print(f"Warning: '{clade}' was not found in the local NCBI taxonomy database.")
        return set()

    clade_taxid = translation[clade][0]

    try:
        descendants = ncbi.get_descendant_taxa(
            clade_taxid,
            intermediate_nodes=True,
        )
    except Exception as error:
        print(f"Warning: could not retrieve descendants of '{clade}': {error}")
        descendants = []

    return {clade_taxid, *descendants}


def make_count_table(factor_taxids, clades):
    ncbi = NCBITaxa()
    clade_taxids = {
        clade: get_clade_taxids(ncbi, clade)
        for clade in clades
    }

    counts = {
        clade: {
            factor: len(factor_taxids[factor].intersection(clade_taxids[clade]))
            for factor in FACTOR_ORDER
        }
        for clade in clades
    }

    count_table = pd.DataFrame.from_dict(counts, orient="index")
    count_table = count_table.loc[clades, FACTOR_ORDER]
    count_table.index.name = "Clade"
    return count_table


def abundance_color(count, normalizer):
    if count <= 0:
        return "white"
    return COLOR_MAP(normalizer(float(count)))


def plot_count_table(count_table, output, dpi, color_max, figure_size):
    if color_max <= 1:
        raise ValueError("--color-max must be greater than 1.")

    normalizer = mpl.colors.Normalize(vmin=1, vmax=color_max, clip=True)
    n_rows, n_columns = count_table.shape
    x_spacing = 2.25
    y_spacing = 1.25
    cell_width = 0.92
    cell_height = 0.82

    figure, axis = plt.subplots(figsize=figure_size)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    for row_index, clade in enumerate(count_table.index):
        y_position = row_index * y_spacing

        for column_index, factor in enumerate(count_table.columns):
            count = count_table.loc[clade, factor]
            x_position = column_index * x_spacing

            axis.add_patch(
                Rectangle(
                    (
                        x_position - cell_width / 2,
                        y_position - cell_height / 2,
                    ),
                    width=cell_width,
                    height=cell_height,
                    facecolor=abundance_color(count, normalizer),
                    edgecolor="#b7b7b7" if count <= 0 else "#8a5a16",
                    linewidth=0.8,
                    zorder=2,
                )
            )

    x_positions = np.arange(n_columns) * x_spacing
    y_positions = np.arange(n_rows) * y_spacing

    axis.set_xticks(x_positions)
    axis.set_xticklabels(count_table.columns, fontsize=11, ha="center")
    axis.xaxis.tick_top()
    axis.tick_params(
        axis="x",
        top=False,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        pad=10,
    )

    axis.set_yticks(y_positions)
    axis.set_yticklabels(count_table.index, fontsize=11)
    axis.tick_params(axis="y", left=False, right=False, pad=8)
    axis.set_xlim(-0.85, (n_columns - 1) * x_spacing + 0.85)
    axis.set_ylim((n_rows - 0.5) * y_spacing, -0.65 * y_spacing)
    axis.set_aspect("auto")
    axis.grid(False)

    for spine in axis.spines.values():
        spine.set_visible(False)

    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main():
    args = get_args()
    factor_taxids = load_factor_taxids(args.input_dir)

    print("Unique input taxids:")
    for factor in FACTOR_ORDER:
        print(f"  {factor}: {len(factor_taxids[factor])}")

    count_table = make_count_table(factor_taxids, args.clades)

    args.counts_output.parent.mkdir(parents=True, exist_ok=True)
    count_table.to_csv(args.counts_output)
    plot_count_table(
        count_table,
        args.output,
        args.dpi,
        args.color_max,
        (args.fig_width, args.fig_height),
    )

    print("\nFungal processing-factor counts:")
    print(count_table.to_string())
    print(f"\nSaved {args.counts_output}")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

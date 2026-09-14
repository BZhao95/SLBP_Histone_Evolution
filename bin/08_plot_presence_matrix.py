#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
from matplotlib.patches import Rectangle


CLADES = [
    "Bilateria",
    "Cnidaria",
    "Placozoa",
    "Ctenophora",
    "Porifera",
    "Choanoflagellata",
    "Filasterea",
    "Ichthyosporea",
    "Fungi",
    "Aphelida",
    "Amoebozoa",
    "Streptophyta",
    "Chlorophyta",
    "Prasinodermophyta",
    "Stramenopiles",
    "Alveolata",
    "Rhizaria",
    "Haptista",
    "Metamonada",
    "Discoba",
    "Rhodophyta",
    "Cryptophyceae",
    "Glaucocystophyceae",
    "Eukaryota incertae sedis",
    "Apusozoa",
    "Hemimastigophora",
    "Provora",
]

FACTOR_FILES = {
    "SLBP": Path("slbp/slbp_taxids.txt"),
    "LSM10": Path("lsm10/lsm10_taxids.txt"),
    "LSM11": Path("lsm11/lsm11_taxids.txt"),
    "U7": Path("U7/u7_taxids.txt"),
    "FLASH": Path("FLASH/flash_taxids.txt"),
    "SLIP1": Path("SLIP1/slip1_taxids.txt"),
    "eIF4G": Path("eIF4G/eif4g_taxids.txt"),
}

COLOR_MAP = mpl.colors.LinearSegmentedColormap.from_list(
    "factor_abundance",
    ["#dbeaf6", "#8dc0ef", "#5397ce", "#246ba3", "#054b89", "#073763"],
)
COLOR_NORMALIZER = mpl.colors.Normalize(vmin=1, vmax=100, clip=True)


def get_args():
    parser = argparse.ArgumentParser(
        description="Count histone-processing factors by clade and plot the results."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing the factor subdirectories and taxid files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("protein_presence_abundance_wide.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "-c",
        "--counts-output",
        type=Path,
        help="Output count table. By default it is saved inside the input directory.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    return parser.parse_args()


def load_taxids(file_path):
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


def get_clade_taxids(ncbi, clade):
    translation = ncbi.get_name_translator([clade])

    if clade not in translation:
        print(f"Warning: '{clade}' was not found in NCBI Taxonomy.")
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


def make_count_table(input_dir):
    factor_taxids = {
        factor: load_taxids(input_dir / relative_path)
        for factor, relative_path in FACTOR_FILES.items()
    }

    ncbi = NCBITaxa()
    clade_taxids = {
        clade: get_clade_taxids(ncbi, clade)
        for clade in CLADES
    }

    counts = {
        factor: {
            clade: len(taxids.intersection(clade_taxids[clade]))
            for clade in CLADES
        }
        for factor, taxids in factor_taxids.items()
    }

    count_table = pd.DataFrame.from_dict(counts, orient="index")
    count_table = count_table.loc[list(FACTOR_FILES), CLADES]
    count_table.index.name = "Factor"
    return count_table


def abundance_color(count):
    if count <= 0:
        return "white"
    return COLOR_MAP(COLOR_NORMALIZER(min(float(count), 100)))


def row_background(clade):
    if clade == "Bilateria":
        return "#F8F2FC"
    if clade == "Fungi":
        return "#FBF6E8"
    if clade in {"Chlorophyta", "Streptophyta", "Prasinodermophyta"}:
        return "#F1FBE9"
    if clade in {"Alveolata", "Stramenopiles", "Rhizaria"}:
        return "#E8EEF6"
    return "white"


def plot_count_table(count_table, output_path, dpi):
    data = count_table.T
    n_rows, n_columns = data.shape

    x_spacing = 4
    y_spacing = 1.5
    square_width = 2
    square_height = 1

    figure, axis = plt.subplots(figsize=(44, 50))

    for row_index, clade in enumerate(data.index):
        y_position = row_index * y_spacing

        axis.add_patch(
            Rectangle(
                (-square_width, y_position - y_spacing / 2),
                width=(n_columns - 1) * x_spacing + 2 * square_width,
                height=y_spacing,
                facecolor=row_background(clade),
                edgecolor="none",
                zorder=0,
            )
        )

        for column_index, factor in enumerate(data.columns):
            count = data.loc[clade, factor]
            x_position = column_index * x_spacing

            axis.add_patch(
                Rectangle(
                    (
                        x_position - square_width / 2,
                        y_position - square_height / 2,
                    ),
                    width=square_width,
                    height=square_height,
                    facecolor=abundance_color(count),
                    edgecolor="#b8b8b8" if count <= 0 else "#4d4d4d",
                    linewidth=1.0,
                    zorder=2,
                )
            )

    x_positions = np.arange(n_columns) * x_spacing
    y_positions = np.arange(n_rows) * y_spacing

    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        data.columns,
        fontsize=16,
        fontweight="bold",
        rotation=45,
        ha="left",
    )
    axis.xaxis.tick_top()
    axis.tick_params(
        axis="x",
        top=False,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        pad=12,
    )

    axis.set_yticks(y_positions)
    axis.set_yticklabels(data.index, fontsize=14)
    axis.tick_params(axis="y", left=False, right=False, pad=10)
    axis.set_xlim(-square_width, (n_columns - 1) * x_spacing + square_width)
    axis.set_ylim((n_rows - 0.5) * y_spacing, -0.5 * y_spacing)
    axis.set_aspect("auto")
    axis.grid(False)

    for spine in axis.spines.values():
        spine.set_visible(False)

    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main():
    args = get_args()
    counts_output = args.counts_output or args.input_dir / "protein_counts_by_clade.csv"

    count_table = make_count_table(args.input_dir)

    counts_output.parent.mkdir(parents=True, exist_ok=True)
    count_table.to_csv(counts_output)
    plot_count_table(count_table, args.output, args.dpi)

    print(count_table.to_string())
    print(f"Saved {counts_output}")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

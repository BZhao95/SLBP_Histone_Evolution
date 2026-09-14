#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba


DEFAULT_TEMPLATE_SL = "CCAAAGGCTCTTTTCAGAGCCACCCA"

HISTONE_MAP = {
    "H2A": ["H2A", "H2A; Provisional", "HFD_H2A"],
    "H2B": ["H2B", "H2B; Provisional", "HFD_H2B"],
    "H3": ["H3", "H3; Provisional", "HFD_H3"],
    "H4": ["H4", "H4; Provisional", "HFD_H4"],
}

PRESETS = {
    "fungi": {
        "clades": [
            "Microsporidia",
            "Cryptomycota",
            "Chytridiomycota",
            "Blastocladiomycota",
            "Zoopagomycota",
            "Mucoromycota",
            "Ascomycota",
            "Basidiomycota",
        ],
        "low_color": "#FFF2CC",
        "high_color": "#7A4A0B",
        "low_alpha": 0.45,
        "high_alpha": 0.90,
        "negative_color": "white",
        "solid_width": 2.4,
        "small_font": 15,
    },
    "metazoan": {
        "clades": [
            "Porifera",
            "Ctenophora",
            "Cnidaria",
            "Placozoa",
            "Bilateria",
        ],
        "low_color": "#E6D5F7",
        "high_color": "#4B006E",
        "low_alpha": 0.35,
        "high_alpha": 0.80,
        "negative_color": "#BFBFBF",
        "solid_width": 2.2,
        "small_font": 18,
    },
    "protozoa": {
        "clades": [
            "Ichthyosporea",
            "Discoba",
            "Metamonada",
            "Amoebozoa",
            "Rhizaria",
            "Stramenopiles",
            "Alveolata",
        ],
        "low_color": "#DCEAF7",
        "high_color": "#1F4E79",
        "low_alpha": 0.40,
        "high_alpha": 0.90,
        "negative_color": "#AFAFAF",
        "solid_width": 2.2,
        "small_font": 15,
    },
    "plant": {
        "clades": [
            "Prasinodermophyta",
            "Chlorophyta",
            "Klebsormidiophyceae",
            "Bryophyta",
            "Acrogymnospermae",
            "eudicotyledons",
            "Liliopsida",
        ],
        "low_color": "#EDF5EE",
        "high_color": "#1E6B3A",
        "low_alpha": 0.50,
        "high_alpha": 1.00,
        "negative_color": "white",
        "solid_width": 2.4,
        "small_font": 15,
    },
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Plot histone stem-loop features across taxonomic clades."
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input TSV file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("histone_sl_bubbles.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--preset",
        choices=PRESETS,
        required=True,
        help="Clade order and colour palette.",
    )
    parser.add_argument(
        "--species-column",
        help="Species column. If omitted, Species_Name or Species is detected.",
    )
    parser.add_argument(
        "--clades",
        nargs="+",
        help="Optional clade order overriding the selected preset.",
    )
    parser.add_argument(
        "--template-sl",
        default=DEFAULT_TEMPLATE_SL,
        help="Template stem-loop sequence used for conservation scoring.",
    )
    parser.add_argument("--core-start", type=int, default=5)
    parser.add_argument("--core-end", type=int, default=21)
    parser.add_argument("--low-color", help="Override the low-conservation colour.")
    parser.add_argument("--high-color", help="Override the high-conservation colour.")
    parser.add_argument("--negative-color", help="Override the SLBP-negative colour.")
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--summary-output",
        type=Path,
        help="Optional TSV containing the values represented by each bubble.",
    )
    parser.add_argument(
        "--network-output",
        type=Path,
        help="Optional top stem-loop network figure.",
    )
    parser.add_argument(
        "--network-summary-output",
        type=Path,
        help="Optional TSV describing the nodes in the stem-loop network.",
    )
    parser.add_argument("--network-top-n", type=int, default=50)
    parser.add_argument("--highlight-clade", default="Bilateria")
    parser.add_argument("--highlight-color", default="purple")
    parser.add_argument("--other-color", default="#B6D7E4")
    parser.add_argument("--network-dpi", type=int, default=300)
    parser.add_argument(
        "--correlation-output",
        type=Path,
        help="Optional conservation-versus-distance figure.",
    )
    parser.add_argument("--distance-min", type=float, default=20)
    parser.add_argument("--distance-max", type=float, default=200)
    parser.add_argument("--distance-bin-size", type=float, default=10)
    parser.add_argument("--correlation-dpi", type=int, default=300)
    return parser.parse_args()


def find_species_column(df, requested_column=None):
    if requested_column:
        if requested_column not in df.columns:
            raise ValueError(f"Species column not found: {requested_column}")
        return requested_column

    for column in ("Species_Name", "Species"):
        if column in df.columns:
            return column

    raise ValueError("Could not find a Species_Name or Species column.")


def conservation_score(sequence, template_core, template_length, core_start, core_end):
    if not isinstance(sequence, str) or len(sequence) != template_length:
        return np.nan

    sequence_core = sequence[core_start:core_end].upper()
    return sum(
        sequence_base == template_base
        for sequence_base, template_base in zip(sequence_core, template_core)
    )


def prepare_data(df, clades, species_column, template_sl, core_start, core_end):
    required_columns = {
        "Clade",
        species_column,
        "CDD_Shortname",
        "SL",
        "SLBP",
        "HDE",
        "SL_Position",
        "SL_Sequence",
        "Protein_Accession",
    }
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    template_sl = template_sl.upper()
    if not 0 <= core_start < core_end <= len(template_sl):
        raise ValueError("The core coordinates fall outside the template sequence.")

    template_core = template_sl[core_start:core_end]
    df = df.copy()
    df["Conservation_to_template"] = df["SL_Sequence"].apply(
        lambda sequence: conservation_score(
            sequence,
            template_core,
            len(template_sl),
            core_start,
            core_end,
        )
    )

    observed_clades = set(df["Clade"].dropna())
    missing_clades = [clade for clade in clades if clade not in observed_clades]

    if missing_clades:
        print(f"Warning: missing clades: {', '.join(missing_clades)}")

    species_per_clade = df.groupby("Clade")[species_column].nunique().to_dict()
    records = []

    for clade in clades:
        clade_df = df[df["Clade"] == clade]
        n_species = species_per_clade.get(clade, 1)

        for histone, names in HISTONE_MAP.items():
            pattern = "|".join(re.escape(name) for name in names)
            histone_df = clade_df[
                clade_df["CDD_Shortname"].str.contains(
                    pattern,
                    case=False,
                    na=False,
                    regex=True,
                )
            ]
            sl_hits = histone_df[
                histone_df["SL"].astype(str).str.upper().eq("Y")
            ]

            n_accessions = sl_hits["Protein_Accession"].nunique()
            normalized_count = n_accessions / n_species if n_species else 0.0
            positions = pd.to_numeric(sl_hits["SL_Position"], errors="coerce")
            conservation = pd.to_numeric(
                sl_hits["Conservation_to_template"],
                errors="coerce",
            )

            records.append(
                {
                    "Clade": clade,
                    "Histone": histone,
                    "N_species": n_species,
                    "N_SL_accessions": n_accessions,
                    "Normalized_SL_accessions": normalized_count,
                    "Median_SL_position": positions.median(),
                    "Mean_core_conservation": conservation.mean(),
                    "SLBP_present": histone_df["SLBP"]
                    .astype(str)
                    .str.upper()
                    .eq("Y")
                    .any(),
                    "HDE_present": histone_df["HDE"]
                    .astype(str)
                    .str.upper()
                    .eq("Y")
                    .any(),
                }
            )

    return pd.DataFrame(records), clades, df


def scale_size(value, minimum, maximum, min_size=200, max_size=1800, power=0.5):
    if value <= 0:
        return min_size * 0.3
    if maximum <= minimum:
        normalized = 0.5
    else:
        normalized = np.clip((value - minimum) / (maximum - minimum), 0, 1)
    return min_size + normalized**power * (max_size - min_size)


def interpolate_color(value, minimum, maximum, low_color, high_color):
    if pd.isna(value):
        return tuple(low_color)
    if maximum <= minimum:
        normalized = 0.5
    else:
        normalized = np.clip((value - minimum) / (maximum - minimum), 0, 1)
    return tuple(low_color + normalized * (high_color - low_color))


def plot_bubbles(summary, clades, preset, output, dpi):
    histones = list(HISTONE_MAP)
    positive_counts = summary.loc[
        summary["Normalized_SL_accessions"] > 0,
        "Normalized_SL_accessions",
    ]
    count_min = positive_counts.min() if not positive_counts.empty else 0.001
    count_max = positive_counts.max() if not positive_counts.empty else 1.0

    conservation_values = summary["Mean_core_conservation"].fillna(0)
    conservation_min = conservation_values.min()
    conservation_max = conservation_values.max()

    low_color = np.array(to_rgba(preset["low_color"], alpha=preset["low_alpha"]))
    high_color = np.array(
        to_rgba(preset["high_color"], alpha=preset["high_alpha"])
    )

    figure, axis = plt.subplots(
        figsize=(3.0 + len(histones) * 2.3, 2.0 + len(clades) * 1.2)
    )

    for y_position, clade in enumerate(clades):
        for x_position, histone in enumerate(histones):
            row = summary[
                summary["Clade"].eq(clade) & summary["Histone"].eq(histone)
            ].iloc[0]

            size = scale_size(
                row["Normalized_SL_accessions"],
                count_min,
                count_max,
            )
            color = (
                interpolate_color(
                    row["Mean_core_conservation"],
                    conservation_min,
                    conservation_max,
                    low_color,
                    high_color,
                )
                if row["SLBP_present"]
                else preset["negative_color"]
            )

            circle = plt.Circle(
                (x_position, y_position),
                radius=np.sqrt(size) / 90,
                facecolor=color,
                edgecolor="black",
                linewidth=preset["solid_width"] if row["HDE_present"] else 1.8,
                linestyle="solid" if row["HDE_present"] else (0, (4, 2)),
            )
            axis.add_patch(circle)

            if pd.notna(row["Median_SL_position"]):
                mean_rgb = np.mean(to_rgba(color)[:3])
                axis.text(
                    x_position,
                    y_position,
                    f"{int(round(row['Median_SL_position']))}",
                    ha="center",
                    va="center",
                    fontsize=20 if size > 700 else preset["small_font"],
                    fontweight="bold",
                    color="white" if mean_rgb < 0.45 else "black",
                )

    axis.set_xticks(np.arange(len(histones)))
    axis.set_xticklabels(histones, fontsize=14, fontweight="bold")
    axis.set_yticks(np.arange(len(clades)))
    axis.set_yticklabels(clades, fontsize=12)
    axis.invert_yaxis()
    axis.set_xlim(-0.6, len(histones) - 0.4)
    axis.set_ylim(-0.6, len(clades) - 0.4)
    axis.set_aspect("equal")
    axis.axis("off")

    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def format_stem_loop(sequence):
    return f"{sequence[:6]} {sequence[6:10]} {sequence[10:]}"


def get_hairpin_energy(
    sequence,
    temperature=27,
    method="nupack mfe",
    dot_paren_structure=False,
    model=None,
):
    try:
        import nupack
    except ImportError as error:
        raise ImportError(
            "Hairpin-energy calculation requires the nupack package."
        ) from error

    sequence = sequence.replace(" ", "")
    if model is None:
        model = nupack.Model(material="rna", celsius=temperature)

    if method == "nupack mfe" and not dot_paren_structure:
        return nupack.mfe(sequence, model=model)[0].energy
    if method == "nupack mfe" and dot_paren_structure:
        return nupack.mfe(sequence, model=model)[0].structure
    if method == "nupack pfunc":
        return nupack.pfunc(sequence, model=model)[1]

    raise ValueError(
        f"Method '{method}' is unavailable; use 'nupack mfe' or 'nupack pfunc'."
    )


def plot_stem_loop_network(
    df,
    output,
    summary_output,
    top_n,
    highlight_clade,
    highlight_color,
    other_color,
    core_start,
    core_end,
    dpi,
):
    try:
        import networkx as nx
    except ImportError as error:
        raise ImportError(
            "Network plotting requires networkx. Install it with: pip install networkx"
        ) from error

    try:
        import stem_loop_distance as sld
    except ImportError as error:
        raise ImportError(
            "Network plotting requires stem_loop_distance.py to be importable."
        ) from error

    if top_n < 1:
        raise ValueError("--network-top-n must be at least 1.")
    if core_end - core_start != 16:
        raise ValueError("The stem-loop network requires a 16-nt core region.")

    network_df = df[["SL_Sequence", "Clade"]].dropna(subset=["SL_Sequence"]).copy()
    network_df["SL_Sequence"] = (
        network_df["SL_Sequence"]
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace("T", "U", regex=False)
    )
    network_df = network_df[network_df["SL_Sequence"].ne("")]
    network_df = network_df[
        network_df["SL_Sequence"].apply(
            lambda sequence: all(base in "ACGU" for base in sequence)
        )
    ].copy()
    network_df["StemLoopCore"] = network_df["SL_Sequence"].str.slice(
        core_start,
        core_end,
    )

    if network_df.empty:
        raise ValueError("No valid stem-loop sequences are available for the network.")

    frequencies = network_df["StemLoopCore"].value_counts()
    top_sequences = frequencies.head(top_n).index.tolist()
    all_sequences = frequencies.index.tolist()

    full_graph = sld.build_graph_from_stem_loops(
        [format_stem_loop(sequence) for sequence in all_sequences],
        distance_func=sld.loop_edit_distance,
    )
    full_graph = nx.relabel_nodes(
        full_graph,
        {
            index: f"{index}:{get_hairpin_energy(sequence):.2f}"
            for index, sequence in enumerate(all_sequences)
        },
    )
    print(f"Full network nodes: {full_graph.number_of_nodes()}")
    print(f"Full network edges: {full_graph.number_of_edges()}")

    formatted_sequences = [format_stem_loop(sequence) for sequence in top_sequences]
    graph = sld.build_graph_from_stem_loops(
        formatted_sequences,
        distance_func=sld.loop_edit_distance,
    )
    graph = nx.relabel_nodes(
        graph,
        {
            index: f"{index}:{get_hairpin_energy(sequence):.2f}"
            for index, sequence in enumerate(top_sequences)
        },
    )

    sequence_clades = (
        network_df.groupby("StemLoopCore")["Clade"]
        .apply(lambda values: set(values.dropna().unique()))
        .to_dict()
    )

    node_colors = []
    node_rows = []

    for node in graph.nodes:
        node_text = str(node).split(":", maxsplit=1)[0]
        if not node_text.isdigit() or int(node_text) >= len(top_sequences):
            raise ValueError(f"Unexpected network node identifier: {node}")

        sequence_index = int(node_text)
        sequence = top_sequences[sequence_index]
        clades = sequence_clades.get(sequence, set())
        is_exclusive = clades == {highlight_clade}
        energy = get_hairpin_energy(sequence)

        node_colors.append(highlight_color if is_exclusive else other_color)
        node_rows.append(
            {
                "Rank": sequence_index + 1,
                "StemLoopCore": sequence,
                "Frequency": int(frequencies[sequence]),
                "Hairpin_energy": energy,
                "Clades": ";".join(sorted(clades)),
                f"Exclusive_to_{highlight_clade}": is_exclusive,
                "Network_degree": graph.degree[node],
            }
        )

    figure = plt.figure(figsize=(10, 8))
    nx.draw_kamada_kawai(
        graph,
        with_labels=False,
        node_size=800,
        node_color=node_colors,
        edge_color="gray",
    )

    x_limits = plt.xlim()
    y_limits = plt.ylim()
    plt.text(
        (x_limits[0] + x_limits[1]) / 2,
        y_limits[0] - (y_limits[1] - y_limits[0]) * 0.08,
        "Average Network Degree: 25",
        fontsize=35,
        ha="center",
        va="top",
    )
    plt.scatter([], [], color=highlight_color, s=120, label=highlight_clade)
    plt.scatter([], [], color=other_color, s=120, label="Others/Shared")
    plt.legend(frameon=False, fontsize=30, loc="upper left")

    output.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)

    if summary_output:
        summary_output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(node_rows).sort_values("Rank").to_csv(
            summary_output,
            sep="\t",
            index=False,
        )

    return graph.number_of_nodes(), graph.number_of_edges()


def plot_conservation_distance(
    df,
    output,
    distance_min,
    distance_max,
    bin_size,
    dpi,
):
    try:
        from scipy.stats import pearsonr, spearmanr
    except ImportError as error:
        raise ImportError("Correlation plotting requires scipy.") from error

    if distance_min >= distance_max:
        raise ValueError("--distance-min must be smaller than --distance-max.")
    if bin_size <= 0:
        raise ValueError("--distance-bin-size must be positive.")

    correlation_df = df[["SL_Position", "Conservation_to_template"]].copy()
    correlation_df["SL_Position"] = pd.to_numeric(
        correlation_df["SL_Position"],
        errors="coerce",
    )
    correlation_df = correlation_df.dropna()
    correlation_df = correlation_df[
        correlation_df["SL_Position"].between(distance_min, distance_max)
    ].copy()

    if correlation_df.empty:
        raise ValueError("No valid observations remain for the correlation plot.")

    correlation_df["Distance_bin"] = (
        np.floor(correlation_df["SL_Position"] / bin_size) * bin_size
    )
    binned = (
        correlation_df.groupby("Distance_bin", as_index=False)[
            "Conservation_to_template"
        ]
        .mean()
        .rename(columns={"Conservation_to_template": "Mean_conservation"})
    )

    if len(binned) < 2:
        raise ValueError("At least two populated distance bins are required.")

    pearson_r, pearson_p = pearsonr(
        binned["Distance_bin"],
        binned["Mean_conservation"],
    )
    spearman_r, spearman_p = spearmanr(
        binned["Distance_bin"],
        binned["Mean_conservation"],
    )

    jitter = (np.random.rand(len(correlation_df)) - 0.5) * 0.2
    slope, intercept = np.polyfit(
        binned["Distance_bin"],
        binned["Mean_conservation"],
        1,
    )
    x_fit = np.linspace(
        binned["Distance_bin"].min(),
        binned["Distance_bin"].max(),
        100,
    )

    with plt.rc_context(
        {
            "font.size": 20,
            "axes.labelsize": 20,
            "axes.titlesize": 20,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
            "legend.fontsize": 18,
        }
    ):
        figure, axis = plt.subplots(figsize=(10, 8))
        axis.scatter(
            correlation_df["SL_Position"],
            correlation_df["Conservation_to_template"] + jitter,
            s=30,
            alpha=0.5,
            color="#696969",
            edgecolors="none",
        )
        axis.scatter(
            binned["Distance_bin"],
            binned["Mean_conservation"],
            s=80,
            color="#4169E1",
            edgecolors="black",
            linewidths=0.5,
        )
        axis.plot(
            x_fit,
            slope * x_fit + intercept,
            color="red",
            linestyle="--",
            linewidth=2,
        )
        axis.set_xlabel("SL distance to stop codon (nt)", labelpad=10)
        axis.set_ylabel("Conservation to reference SL", labelpad=10)
        axis.set_ylim(-1, 17)
        axis.text(
            0.95,
            0.95,
            f"Pearson r = {pearson_r:.3f}, p = {pearson_p:.1e}\n"
            f"Spearman ρ = {spearman_r:.3f}, p = {spearman_p:.1e}",
            transform=axis.transAxes,
            fontsize=18,
            ha="right",
            va="top",
            bbox={"facecolor": "white", "alpha": 0.6, "edgecolor": "none"},
        )
        figure.tight_layout()
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=dpi, bbox_inches="tight")
        plt.close(figure)

    return pearson_r, pearson_p, spearman_r, spearman_p, len(correlation_df), len(binned)


def main():
    args = get_args()
    preset = PRESETS[args.preset].copy()
    preset["low_color"] = args.low_color or preset["low_color"]
    preset["high_color"] = args.high_color or preset["high_color"]
    preset["negative_color"] = args.negative_color or preset["negative_color"]

    df = pd.read_csv(args.input, sep="\t", dtype={"SL_Sequence": str})
    species_column = find_species_column(df, args.species_column)
    clades = args.clades or preset["clades"]

    summary, available_clades, scored_df = prepare_data(
        df,
        clades,
        species_column,
        args.template_sl,
        args.core_start,
        args.core_end,
    )
    plot_bubbles(summary, available_clades, preset, args.output, args.dpi)

    if args.summary_output:
        args.summary_output.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary_output, sep="\t", index=False)
        print(f"Saved {args.summary_output}")

    print(f"Saved {args.output}")

    if args.network_summary_output and not args.network_output:
        raise ValueError("--network-summary-output requires --network-output.")

    if args.network_output:
        n_nodes, n_edges = plot_stem_loop_network(
            scored_df,
            args.network_output,
            args.network_summary_output,
            args.network_top_n,
            args.highlight_clade,
            args.highlight_color,
            args.other_color,
            args.core_start,
            args.core_end,
            args.network_dpi,
        )
        print(f"Top network: {n_nodes} nodes, {n_edges} edges")
        if args.network_summary_output:
            print(f"Saved {args.network_summary_output}")
        print(f"Saved {args.network_output}")

    if args.correlation_output:
        pearson_r, pearson_p, spearman_r, spearman_p, n_rows, n_bins = (
            plot_conservation_distance(
                scored_df,
                args.correlation_output,
                args.distance_min,
                args.distance_max,
                args.distance_bin_size,
                args.correlation_dpi,
            )
        )
        print(f"Correlation: {n_rows} observations in {n_bins} bins")
        print(f"Pearson r = {pearson_r:.3f}, p = {pearson_p:.1e}")
        print(f"Spearman ρ = {spearman_r:.3f}, p = {spearman_p:.1e}")
        print(f"Saved {args.correlation_output}")


if __name__ == "__main__":
    main()

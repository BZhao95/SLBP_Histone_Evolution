#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


COL_SPEC = "Species"
COL_CLADE = "Clade"
COL_SLBP = "SLBP"
COL_SL = "SL"
COL_PAS = "PAS"

CLADES = [
    "Cryptomycota",
    "Blastocladiomycota",
    "Chytridiomycota",
    "Zoopagomycota",
    "Microsporidia",
    "Mucoromycota",
]

STATE_ORDER = ["SLBP–SL", "SLBP–SL–PAS", "SLBP–PAS", "PAS-only"]
STATE_COLORS = {
    "SLBP–SL": "#009E73",
    "SLBP–SL–PAS": "#0072B2",
    "SLBP–PAS": "#E69F00",
    "PAS-only": "#CC79A7",
}

RIDGE_DRAW_ORDER = STATE_ORDER.copy()
SPECIES_EXHIBIT_THRESHOLD = 0.0
SIGMA_DDOF = 0

RIDGE_TAIL_CUTOFF = 0.005
OUTLINE_VISIBILITY_TOLERANCE = 1e-4
IDENTICAL_LAYER_LINEWIDTH = 0.55
NORMAL_FILL_ALPHA = 0.96
OVERLAPPING_EQUAL_HEIGHT_ALPHA = 0.72
EQUAL_HEIGHT_TOLERANCE = 1e-8
STRONG_OVERLAP_SIGMA_FACTOR = 1.0
MARKER_SIZE = 5.5
MARKER_EDGE_WIDTH = 0.85


def output_path(input_path, suffix):
    return input_path.with_name(f"{input_path.stem}{suffix}")


def is_yes(value):
    return str(value).strip().upper() == "Y"


def has_pas(value):
    value = str(value).strip()
    return bool(value) and value.upper() != "NA"


def classify_row(row):
    slbp, sl, pas = row["_SLBP"], row["_SL"], row["_PAS"]

    if slbp and sl and not pas:
        return "SLBP–SL"
    if slbp and sl and pas:
        return "SLBP–SL–PAS"
    if slbp and not sl and pas:
        return "SLBP–PAS"
    if not slbp and not sl and pas:
        return "PAS-only"
    return "Other"


def dominant_state(row):
    values = [row[state] for state in STATE_ORDER]

    if sum(values) == 0:
        return "Other"

    maximum = max(values)
    ties = [state for state, value in zip(STATE_ORDER, values) if value == maximum]
    return "SLBP–SL–PAS" if "SLBP–SL–PAS" in ties else ties[0]


def safe_logit(values):
    values = np.clip(values, 1e-6, 1 - 1e-6)
    return np.log(values / (1 - values))


def bootstrap_mean_lor(values_sl, values_pas, rng, n_boot=3000):
    lor_per_species = safe_logit(values_sl) - safe_logit(values_pas)

    if len(lor_per_species) == 0:
        return np.nan, np.nan, np.nan

    indices = rng.integers(
        0,
        len(lor_per_species),
        size=(n_boot, len(lor_per_species)),
    )
    bootstrap_means = lor_per_species[indices].mean(axis=1)
    low, high = np.percentile(bootstrap_means, [2.5, 97.5])

    return float(bootstrap_means.mean()), float(low), float(high)


def prepare_data(input_path):
    df = pd.read_csv(input_path, sep="\t", dtype=str).fillna("NA")

    required_columns = {COL_SPEC, COL_CLADE, COL_SLBP, COL_SL, COL_PAS}
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    df["_SLBP"] = df[COL_SLBP].apply(is_yes)
    df["_SL"] = df[COL_SL].apply(is_yes)
    df["_PAS"] = df[COL_PAS].apply(has_pas)
    df["State"] = df.apply(classify_row, axis=1)

    df_focus = df[df["State"].isin(STATE_ORDER)].copy()
    if df_focus.empty:
        raise ValueError("No loci matched the four processing states.")

    grouped_states = df_focus.groupby([COL_SPEC, COL_CLADE])["State"]
    species_comp = (
        grouped_states.value_counts(normalize=True)
        .unstack(fill_value=0)
        .reindex(columns=STATE_ORDER, fill_value=0)
    )
    species_comp["N_loci"] = grouped_states.size()
    species_comp = species_comp.reset_index().rename(
        columns={COL_SPEC: "Species", COL_CLADE: "Clade"}
    )

    clade_mean = (
        species_comp.groupby("Clade")[STATE_ORDER]
        .mean()
        .reindex(CLADES)
        .fillna(0.0)
    )

    species_comp["Dominant_State"] = species_comp.apply(dominant_state, axis=1)
    dominant_counts = (
        species_comp[species_comp["Dominant_State"].isin(STATE_ORDER)]
        .groupby(["Clade", "Dominant_State"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=CLADES, columns=STATE_ORDER)
        .fillna(0)
    )

    species_comp["_p_SL"] = (
        species_comp["SLBP–SL"] + species_comp["SLBP–SL–PAS"]
    ).astype(float)
    species_comp["_p_PASonly"] = species_comp["PAS-only"].astype(float)

    rng = np.random.default_rng(42)
    log_odds_rows = []

    for clade in CLADES:
        subset = species_comp[species_comp["Clade"] == clade]
        mean, low, high = bootstrap_mean_lor(
            subset["_p_SL"].to_numpy(),
            subset["_p_PASonly"].to_numpy(),
            rng,
        )
        log_odds_rows.append(
            {
                "Clade": clade,
                "mean_LOR": mean,
                "CI_low": low,
                "CI_high": high,
            }
        )

    log_odds = pd.DataFrame(log_odds_rows).set_index("Clade").reindex(CLADES)
    return species_comp, clade_mean, dominant_counts, log_odds


def species_sem(values):
    values = np.asarray(values, dtype=float)

    if len(values) < 2:
        return 0.0, "insufficient_data"

    standard_deviation = float(np.std(values, ddof=SIGMA_DDOF))
    if standard_deviation <= 0:
        return 0.0, "zero_variance"

    return standard_deviation / np.sqrt(len(values)), "ok"


def boundary_capped_sigma(raw_sigma, mean_x):
    boundary_distance = min(mean_x, 1.0 - mean_x)

    if boundary_distance <= 0 or raw_sigma <= boundary_distance:
        return raw_sigma, False

    return boundary_distance, True


def scaled_amplitude(species_proportion, exponent=0.5, maximum=1.0):
    if species_proportion <= 0:
        return 0.0

    amplitude = float(np.clip(species_proportion, 0, 1)) ** exponent
    return min(amplitude, maximum)


def bump_curve(x_values, mean_x, sigma):
    curve = np.exp(-0.5 * ((x_values - mean_x) / sigma) ** 2)
    return curve / curve.max() if curve.max() > 0 else curve


def plot_visible_outline(
    axis,
    x_values,
    curve,
    visible_mask,
    color="black",
    linewidth=0.85,
    zorder=20,
):
    indices = np.flatnonzero(visible_mask)
    if len(indices) == 0:
        return

    breaks = np.where(np.diff(indices) > 1)[0]
    starts = np.r_[indices[0], indices[breaks + 1]]
    ends = np.r_[indices[breaks], indices[-1]]

    for start, end in zip(starts, ends):
        if end - start < 3:
            continue
        axis.plot(
            x_values[start : end + 1],
            curve[start : end + 1],
            color=color,
            linewidth=linewidth,
            zorder=zorder,
            solid_joinstyle="round",
            solid_capstyle="round",
        )


def group_identical_geometry(ridge_curves, draw_order, decimal_places=8):
    grouped = {}

    for state in draw_order:
        ridge = ridge_curves[state]
        key = (
            round(ridge["mean_x"], decimal_places),
            round(ridge["amplitude"], decimal_places),
            round(ridge["sigma"], decimal_places),
        )
        grouped.setdefault(key, []).append(state)

    groups = []
    for states in grouped.values():
        ridge = ridge_curves[states[0]]
        groups.append(
            {
                "states": states,
                "amplitude": ridge["amplitude"],
                "mean_x": ridge["mean_x"],
                "sigma": ridge["sigma"],
            }
        )

    state_position = {state: i for i, state in enumerate(RIDGE_DRAW_ORDER)}
    return sorted(
        groups,
        key=lambda group: (
            -group["amplitude"],
            min(state_position[state] for state in group["states"]),
        ),
    )


def groups_overlap(group_a, group_b):
    equal_height = (
        abs(group_a["amplitude"] - group_b["amplitude"])
        <= EQUAL_HEIGHT_TOLERANCE
    )
    centre_distance = abs(group_a["mean_x"] - group_b["mean_x"])
    overlap_limit = STRONG_OVERLAP_SIGMA_FACTOR * max(
        group_a["sigma"], group_b["sigma"]
    )
    return equal_height and centre_distance < overlap_limit


def plot_panel_a(species_comp, input_path):
    clades_a = list(reversed(CLADES))
    x_values = np.linspace(0, 1, 1601)
    clade_gap = 1.35
    clade_y = {
        clade: (len(clades_a) - 1 - i) * clade_gap
        for i, clade in enumerate(clades_a)
    }

    precomputed = {}
    for clade in clades_a:
        subset = species_comp[species_comp["Clade"] == clade]

        for state in RIDGE_DRAW_ORDER:
            values = subset[state].to_numpy(dtype=float)
            if len(values) == 0:
                continue

            n_exhibiting = int(np.sum(values > SPECIES_EXHIBIT_THRESHOLD))
            species_proportion = n_exhibiting / len(values)
            raw_sem, sem_status = species_sem(values)

            precomputed[(clade, state)] = {
                "mean_x": float(values.mean()),
                "n_total": len(values),
                "n_exhibiting": n_exhibiting,
                "species_proportion": species_proportion,
                "amplitude": scaled_amplitude(species_proportion),
                "raw_sem": raw_sem,
                "sem_status": sem_status,
            }

    measured_sems = [
        values["raw_sem"]
        for values in precomputed.values()
        if values["sem_status"] == "ok" and values["amplitude"] > 0
    ]
    minimum_measured_sem = min(measured_sems) if measured_sems else None

    figure, axis = plt.subplots(figsize=(8.4, 6.4))

    for x_position in [0.2, 0.4, 0.6, 0.8]:
        axis.axvline(
            x_position,
            color="#BDBDBD",
            linestyle=(0, (1.4, 2.4)),
            linewidth=0.9,
            zorder=0,
        )

    overall_top = -np.inf
    bump_rows = []

    for clade in clades_a:
        baseline = clade_y[clade]
        axis.plot([0, 1], [baseline, baseline], color="black", linewidth=0.85)

        ridge_curves = {}
        marker_states = []

        for state in RIDGE_DRAW_ORDER:
            if (clade, state) not in precomputed:
                continue

            values = precomputed[(clade, state)]
            sigma_source = "measured"
            marker_reason = ""
            marker_only = False
            was_capped = False

            if values["sem_status"] == "ok":
                sigma, was_capped = boundary_capped_sigma(
                    values["raw_sem"], values["mean_x"]
                )
            elif (
                values["sem_status"] == "zero_variance"
                and minimum_measured_sem is not None
            ):
                sigma, was_capped = boundary_capped_sigma(
                    minimum_measured_sem, values["mean_x"]
                )
                sigma_source = "borrowed_min_dataset_sem"
            else:
                sigma = np.nan
                sigma_source = "none"
                marker_reason = values["sem_status"]
                marker_only = True

            bump_rows.append(
                {
                    "Clade": clade,
                    "State": state,
                    "Mean_x": values["mean_x"],
                    "N_species_exhibiting": values["n_exhibiting"],
                    "N_species_total": values["n_total"],
                    "Species_proportion": values["species_proportion"],
                    "Bump_height": values["amplitude"],
                    "Sigma_sem_raw": (
                        values["raw_sem"] if values["n_total"] >= 2 else np.nan
                    ),
                    "Sigma_final": sigma,
                    "Sigma_source": sigma_source,
                    "Was_capped": was_capped,
                    "Marker_only": marker_only,
                    "Marker_reason": marker_reason,
                }
            )

            if values["amplitude"] <= 0:
                continue

            if marker_only:
                marker_states.append(
                    {
                        "state": state,
                        "mean_x": values["mean_x"],
                        "amplitude": values["amplitude"],
                        "reason": marker_reason,
                    }
                )
                continue

            height_curve = (
                bump_curve(x_values, values["mean_x"], sigma)
                * values["amplitude"]
            )
            ridge_curves[state] = {
                "height_curve": height_curve,
                "upper": baseline + height_curve,
                "amplitude": values["amplitude"],
                "mean_x": values["mean_x"],
                "sigma": sigma,
            }

        state_position = {state: i for i, state in enumerate(RIDGE_DRAW_ORDER)}
        draw_order = sorted(
            ridge_curves,
            key=lambda state: (
                -ridge_curves[state]["amplitude"],
                state_position[state],
            ),
        )
        geometry_groups = group_identical_geometry(ridge_curves, draw_order)

        group_curves = []
        for draw_index, group in enumerate(geometry_groups):
            ridge = ridge_curves[group["states"][0]]
            group_curves.append(
                {
                    **group,
                    "height_curve": ridge["height_curve"],
                    "upper": ridge["upper"],
                    "fill_mask": ridge["height_curve"] > RIDGE_TAIL_CUTOFF,
                    "draw_index": draw_index,
                }
            )

        for group_index, group in enumerate(group_curves):
            fill_alpha = (
                OVERLAPPING_EQUAL_HEIGHT_ALPHA
                if any(
                    groups_overlap(group, other)
                    for other_index, other in enumerate(group_curves)
                    if other_index != group_index
                )
                else NORMAL_FILL_ALPHA
            )

            if len(group["states"]) == 1:
                axis.fill_between(
                    x_values,
                    baseline,
                    group["upper"],
                    where=group["fill_mask"],
                    interpolate=True,
                    color=STATE_COLORS[group["states"][0]],
                    alpha=fill_alpha,
                    linewidth=0,
                    zorder=2 + group["draw_index"],
                )
            else:
                n_layers = len(group["states"])
                for layer_index, state in enumerate(group["states"]):
                    lower = baseline + group["height_curve"] * layer_index / n_layers
                    upper = (
                        baseline
                        + group["height_curve"] * (layer_index + 1) / n_layers
                    )
                    axis.fill_between(
                        x_values,
                        lower,
                        upper,
                        where=group["fill_mask"],
                        interpolate=True,
                        color=STATE_COLORS[state],
                        alpha=fill_alpha,
                        linewidth=0,
                        zorder=2 + group["draw_index"],
                    )

            if group["fill_mask"].any():
                overall_top = max(
                    overall_top,
                    float(np.max(group["upper"][group["fill_mask"]])),
                )

        for group_index, group in enumerate(group_curves):
            covered = np.zeros_like(x_values, dtype=bool)

            for front_group in group_curves[group_index + 1 :]:
                if groups_overlap(group, front_group):
                    continue
                covered |= (
                    front_group["fill_mask"]
                    & (
                        front_group["upper"]
                        >= group["upper"] - OUTLINE_VISIBILITY_TOLERANCE
                    )
                )

            plot_visible_outline(
                axis,
                x_values,
                group["upper"],
                group["fill_mask"] & ~covered,
                zorder=20 + group_index,
            )

            if len(group["states"]) <= 1:
                continue

            n_layers = len(group["states"])
            for boundary_index in range(1, n_layers):
                internal_curve = (
                    baseline
                    + group["height_curve"] * boundary_index / n_layers
                )
                internal_covered = np.zeros_like(x_values, dtype=bool)

                for front_group in group_curves[group_index + 1 :]:
                    if groups_overlap(group, front_group):
                        continue
                    internal_covered |= (
                        front_group["fill_mask"]
                        & (
                            front_group["upper"]
                            >= internal_curve - OUTLINE_VISIBILITY_TOLERANCE
                        )
                    )

                plot_visible_outline(
                    axis,
                    x_values,
                    internal_curve,
                    group["fill_mask"] & ~internal_covered,
                    linewidth=IDENTICAL_LAYER_LINEWIDTH,
                    zorder=19 + group_index,
                )

        for marker in marker_states:
            y_point = baseline + marker["amplitude"]
            color = STATE_COLORS[marker["state"]]

            if marker["reason"] == "zero_variance":
                axis.plot(
                    [marker["mean_x"], marker["mean_x"]],
                    [baseline, y_point],
                    color=color,
                    linewidth=2.2,
                    solid_capstyle="butt",
                    zorder=30,
                    clip_on=False,
                )
                axis.plot(
                    marker["mean_x"],
                    y_point,
                    marker="D",
                    markersize=MARKER_SIZE,
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markeredgewidth=MARKER_EDGE_WIDTH,
                    zorder=31,
                    clip_on=False,
                )
            else:
                axis.plot(
                    marker["mean_x"],
                    y_point,
                    marker="o",
                    markersize=MARKER_SIZE,
                    markerfacecolor="white",
                    markeredgecolor=color,
                    markeredgewidth=1.4,
                    zorder=30,
                    clip_on=False,
                )

            overall_top = max(overall_top, y_point)

    bump_table = pd.DataFrame(bump_rows)
    bump_table.to_csv(
        output_path(input_path, "_panelA_bump_heights_sem_width.tsv"),
        sep="\t",
        index=False,
    )

    axis.set_xlim(0, 1)
    axis.set_ylim(-0.2, overall_top + 0.35 if np.isfinite(overall_top) else 1.0)
    axis.set_xlabel("Mean per-species fraction of loci", fontsize=13)
    axis.set_xticks(np.linspace(0, 1, 6))
    axis.set_xticklabels([f"{value:.1f}" for value in np.linspace(0, 1, 6)])
    axis.set_yticks([clade_y[clade] for clade in clades_a])
    axis.set_yticklabels([])
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.0)
    axis.tick_params(axis="y", length=0)

    legend_handles = [
        Patch(facecolor=STATE_COLORS[state], edgecolor="none", label=state)
        for state in STATE_ORDER
    ]
    axis.legend(
        handles=legend_handles,
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.13),
        columnspacing=1.4,
        handlelength=1.25,
    )

    figure.tight_layout()
    png_path = output_path(input_path, "_panelA.png")
    pdf_path = output_path(input_path, "_panelA.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def plot_panel_b(dominant_counts, input_path):
    figure, axis = plt.subplots(figsize=(7.2, 4.4))
    x_values = np.arange(len(dominant_counts))
    width = 0.18

    for index, state in enumerate(STATE_ORDER):
        axis.bar(
            x_values + index * width,
            dominant_counts[state].to_numpy(),
            width=width,
            color=STATE_COLORS[state],
            edgecolor="black",
            linewidth=0.85,
            label=state,
            zorder=3,
        )

    axis.set_xticks(x_values + width * (len(STATE_ORDER) - 1) / 2)
    axis.set_xticklabels(dominant_counts.index, rotation=35, ha="right")
    axis.set_ylabel("Species count", fontsize=13)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.7, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.0)
    axis.spines["bottom"].set_linewidth(1.0)
    axis.tick_params(axis="both", labelsize=10)
    axis.legend(
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18),
        columnspacing=1.2,
        handlelength=1.3,
    )

    figure.tight_layout()
    png_path = output_path(input_path, "_panelB_dominant_state_counts.png")
    pdf_path = output_path(input_path, "_panelB_dominant_state_counts.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def plot_panel_c(log_odds, input_path):
    figure, axis = plt.subplots(figsize=(6.0, 4.4))
    x_values = np.arange(len(log_odds))
    mean_lor = log_odds["mean_LOR"].to_numpy()
    error = np.vstack(
        [
            mean_lor - log_odds["CI_low"].to_numpy(),
            log_odds["CI_high"].to_numpy() - mean_lor,
        ]
    )

    axis.errorbar(
        x_values,
        mean_lor,
        yerr=error,
        fmt="o",
        capsize=4,
        color="black",
        ecolor="black",
        elinewidth=1.0,
        markersize=4.5,
        zorder=3,
    )
    axis.axhline(0, color="grey", linestyle="--", linewidth=1.0, zorder=1)
    axis.set_xticks(x_values)
    axis.set_xticklabels(log_odds.index, rotation=35, ha="right")
    axis.set_ylabel("Log-odds (SL-containing vs PAS-only)", fontsize=13)
    axis.grid(False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.0)
    axis.spines["bottom"].set_linewidth(1.0)
    axis.tick_params(axis="both", labelsize=10)

    figure.tight_layout()
    png_path = output_path(input_path, "_panelC_logodds.png")
    pdf_path = output_path(input_path, "_panelC_logodds.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def main():
    parser = argparse.ArgumentParser(
        description="Compare processing states across basal fungal clades."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("data/histone_fungi/fungi_meta.tsv"),
        help="Path to the fungal metadata TSV file.",
    )
    args = parser.parse_args()

    species_comp, clade_mean, dominant_counts, log_odds = prepare_data(args.input)

    species_comp.to_csv(
        output_path(args.input, "_species_state_fractions.tsv"),
        sep="\t",
        index=False,
    )
    clade_mean.to_csv(
        output_path(args.input, "_clade_mean_state_fractions.tsv"),
        sep="\t",
    )

    figure_paths = [
        *plot_panel_a(species_comp, args.input),
        *plot_panel_b(dominant_counts, args.input),
        *plot_panel_c(log_odds, args.input),
    ]

    for path in figure_paths:
        print(f"Saved {path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create publication-ready candidate figures for the Scheme C extension."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SCHEME_C_ROOT", Path(__file__).resolve().parent))
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

ORDER = ["T1", "T2", "T3", "NEURO", "ADHD", "ASD", "ASTHMA"]
DISPLAY = {
    "T1": "T1: affective instability",
    "T2": "T2: help-seeking distress",
    "T3": "T3: major depression",
    "NEURO": "Neuroticism",
    "ADHD": "ADHD",
    "ASD": "Autism spectrum disorder",
    "ASTHMA": "Asthma",
}


def load_pairs() -> pd.DataFrame:
    return pd.read_csv(RESULTS / "Scheme_C_GenomicSEM_pairwise_rg.csv")


def matrix_from_pairs(pairs: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.DataFrame(np.eye(len(ORDER)), index=ORDER, columns=ORDER)
    for row in pairs.itertuples(index=False):
        matrix.loc[row.trait_1, row.trait_2] = row.rg
        matrix.loc[row.trait_2, row.trait_1] = row.rg
    return matrix.rename(index=DISPLAY, columns=DISPLAY)


def heatmap(pairs: pd.DataFrame) -> None:
    matrix = matrix_from_pairs(pairs)
    mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
    fig, ax = plt.subplots(figsize=(8.2, 7.1))
    sns.heatmap(
        matrix,
        mask=mask,
        cmap=sns.diverging_palette(15, 245, as_cmap=True),
        vmin=-1,
        vmax=1,
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 9},
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": r"Genetic correlation ($r_g$)", "shrink": 0.82},
        square=True,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    ax.set_title("Genome-wide genetic correlations across the targeted seven-trait panel", pad=14)
    fig.tight_layout()
    fig.savefig(FIGURES / "Scheme_C_genetic_correlation_heatmap.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "Scheme_C_genetic_correlation_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def tier_profile(pairs: pd.DataFrame) -> None:
    comparators = ["ADHD", "ASD", "NEURO", "ASTHMA"]
    tiers = ["T1", "T2", "T3"]
    colors = {"T1": "#355C7D", "T2": "#C06C84", "T3": "#6C5B7B"}
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.9), sharey=True)
    for ax, comparator in zip(axes, comparators):
        subset_rows = []
        for tier in tiers:
            match = pairs[
                ((pairs.trait_1 == tier) & (pairs.trait_2 == comparator))
                | ((pairs.trait_2 == tier) & (pairs.trait_1 == comparator))
            ]
            if len(match) != 1:
                raise RuntimeError(f"Expected one row for {tier} and {comparator}")
            subset_rows.append(match.iloc[0])
        for y, (tier, row) in enumerate(zip(tiers, subset_rows)):
            ax.errorbar(
                row.rg,
                y,
                xerr=[[row.rg - row.ci_low], [row.ci_high - row.rg]],
                fmt="o",
                color=colors[tier],
                ecolor=colors[tier],
                capsize=3,
                markersize=5,
                linewidth=1.4,
            )
        ax.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_title(DISPLAY[comparator], fontsize=10)
        ax.grid(axis="x", color="#E5E5E5", linewidth=0.7)
        ax.set_axisbelow(True)
    axes[0].set_yticks(range(3), [DISPLAY[tier] for tier in tiers])
    axes[0].invert_yaxis()
    axes[0].set_ylabel("")
    fig.supxlabel(r"Genetic correlation ($r_g$), with 95% confidence interval", y=0.02)
    fig.suptitle("Tier-specific genetic-correlation profiles", y=1.02, fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(FIGURES / "Scheme_C_tier_profile_forest.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "Scheme_C_tier_profile_forest.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    sns.set_theme(style="white", context="paper", font_scale=1.0)
    pairs = load_pairs()
    heatmap(pairs)
    tier_profile(pairs)
    print("Created two PNG and two vector PDF figures")


if __name__ == "__main__":
    main()

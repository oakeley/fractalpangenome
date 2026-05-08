"""
Plotting functions for BlueGenomics
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle
from pathlib import Path
from typing import Union, Dict, Any


def plot_karyoplot(
    genome, annotations: dict, jitter: bool = True, include_mito: bool = False
):
    """
    Generate karyotype plot with integration site annotations.

    Args:
        genome: Genome object with identifier() method and chromosomes attribute
        annotations: Dictionary mapping chromosome names to lists of annotation dicts.
                    Each annotation dict must contain 'pos' key with genomic position.
                    Optional keys: 'offset' (float), 'color' (str)
        jitter: Add random jitter to annotation positions for visualization
        include_mito: Include mitochondrial chromosome in plot

    Example annotations format:
        {
            'chr1': [
                {'pos': 1000000, 'color': 'blue'},
                {'pos': 2000000, 'color': 'red'}
            ],
            'chr2': [...]
        }
    """
    from bluegenomics.config import config

    dataset_dir = config.dataset_directory
    supported_genomes = {
        "hg38": dataset_dir / "Human_hg38_ideogram.tsv",
        "mm10": dataset_dir / "Mouse_mm10_ideogram.tsv",
    }

    # Test for supported genome
    genome_id = genome.identifier()
    assert (
        genome_id in supported_genomes
    ), f"Unsupported genome: {genome_id}, select one of: {list(supported_genomes.keys())}"

    chromosomes = [
        chromosome
        for chromosome in genome.chromosomes
        if "M" not in chromosome or include_mito
    ]

    karyo_dict = {}
    with open(supported_genomes[genome_id]) as karyo_f:
        all_lines = [x.replace(os.linesep, "").split() for x in karyo_f.readlines() if x.strip()]
        lines = [x for x in all_lines if len(x) >= 5 and not x[0].startswith('#')]

        for chromosome in chromosomes:
            karyo_dict[chromosome] = [
                [y[0], int(y[1]), int(y[2]), y[3], y[4]]
                for y in lines if y[0] == chromosome
            ]

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)

    DIM = 1.0

    def get_chromosome_length(chromosome):
        if not karyo_dict[chromosome]:
            return 0
        chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
        chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
        chromosome_length = chromosome_end - chromosome_start
        return chromosome_length

    chromosome_lengths = {
        chromosome: get_chromosome_length(chromosome)
        for chromosome in chromosomes
        if karyo_dict[chromosome]  # Only include chromosomes with band data
    }

    ax.set_xlim([0.0, DIM * (1.1)])
    ax.set_ylim([0.0, DIM])

    def plot_chromosome(chromosome, order):
        chromosome_length = chromosome_lengths[chromosome]
        chromosome_length_1 = max(chromosome_lengths.values())

        x_start = order * DIM * (1 / len(chromosome_lengths))
        x_end = x_start + (DIM * (1 / (len(chromosome_lengths) * 2)))
        y_start = DIM * 0.9 * (chromosome_length / chromosome_length_1)
        y_end = DIM * 0.1

        colors = {
            "gpos100": (0 / 255.0, 0 / 255.0, 0 / 255.0),
            "gpos": (0 / 255.0, 0 / 255.0, 0 / 255.0),
            "gpos75": (130 / 255.0, 130 / 255.0, 130 / 255.0),
            "gpos66": (160 / 255.0, 160 / 255.0, 160 / 255.0),
            "gpos50": (200 / 255.0, 200 / 255.0, 200 / 255.0),
            "gpos33": (210 / 255.0, 210 / 255.0, 210 / 255.0),
            "gpos25": (200 / 255.0, 200 / 255.0, 200 / 255.0),
            "gvar": (220 / 255.0, 220 / 255.0, 220 / 255.0),
            "gneg": (255 / 255.0, 255 / 255.0, 255 / 255.0),
            "acen": (217 / 255.0, 47 / 255.0, 39 / 255.0),
            "stalk": (100 / 255.0, 127 / 255.0, 164 / 255.0),
        }

        for index, piece in enumerate(karyo_dict[chromosome]):
            current_height = piece[2] - piece[1]
            current_height_sc = ((y_end - y_start) / chromosome_length) * current_height
            if index == 0:
                y_previous = y_start

            y_next = y_previous + current_height_sc
            color = colors[piece[4]]

            r = Rectangle(
                (x_start, y_previous),
                x_end - x_start,
                current_height_sc,
                color=color,
            )
            ax.add_patch(r)
            y_previous = y_next

        # Plot semicircles at the beginning and end of the chromosomes
        center_x = x_start + (x_end - x_start) / 2.0
        radius = (x_end - x_start) / 2.0
        theta1 = 0.0
        theta2 = 180.0
        w1 = Wedge(
            (center_x, y_start),
            radius,
            theta1,
            theta2,
            width=0.00001,
            facecolor="white",
            edgecolor="black",
        )
        w2 = Wedge(
            (center_x, y_end),
            radius,
            theta2,
            theta1,
            width=0.00001,
            facecolor="white",
            edgecolor="black",
        )
        ax.add_patch(w1)
        ax.add_patch(w2)
        ax.plot([x_start, x_start], [y_start, y_end], ls="-", color="black")
        ax.plot([x_end, x_end], [y_start, y_end], ls="-", color="black")

        # Plot metadata
        for md in annotations.get(chromosome, []):
            try:
                pos = md["pos"]
            except (KeyError, IndexError):
                err_msg = 'annotations not formatted correctly, annotations must contain "pos" key'
                raise ValueError(err_msg)

            offset = md.get(
                "offset",
                np.random.normal(loc=0.009, scale=0.004) if jitter else 0.015,
            )
            color = md.get("color", "blue")
            ax.plot(
                [x_end + (DIM * 0.003 + offset)],
                [y_start + (y_end - y_start) * (pos / chromosome_length)],
                ".",
                color=color,
            )

        ax.text(center_x, y_end - (DIM * 0.07), chromosome)

    # Only plot chromosomes that have band data
    chromosomes_to_plot = [c for c in chromosomes if c in chromosome_lengths]
    for idx, chromosome in enumerate(chromosomes_to_plot):
        plot_chromosome(chromosome, idx)

    plt.axis("off")
    plt.show()

import edlib
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser
import sys
from collections import defaultdict as d
import matplotlib.pyplot as plt
import matplotlib
import warnings

argparse = ArgumentParser()
argparse.add_argument(
    "-p", "--primers", help="Path to the csv file containing primers", required=True
)
argparse.add_argument(
    "-o", "--output", help="Path to the output csv file", required=True)
argparse.add_argument(
    "-hr",
    "--human_readable",
    help="Avoid annotating heatmap with more than 3 primer pairs",
    required=False,
    default=False,
    action="store_true",
)
argparse.add_argument(
    "-mp",
    "--many_primers",
    help="Avoid creating the heatmap when you have too many primers (>15)",
    required=False,
    default=False,
    action="store_true",
)


args = argparse.parse_args()
prim = args.primers
outf = args.output
humanread = args.human_readable
manyp = args.many_primers


class HumanReadableWarning(UserWarning):
    """Print this warning when there are more than 10 pairs of primers and --human_readable is enabled"""


class ManyPrimersWarning(UserWarning):
    """Print this warning when there are more than 15 pairs of primers and --many_primers is enabled"""


def heatmap(
    data, row_labels, col_labels, ax=None, cbar_kw=None, cbarlabel="", **kwargs
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-90,
             ha="center", rotation_mode="default")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    fontsize=None,
    textcolors=("black", "white"),
    threshold=None,
    **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(
                data[i, j], None), fontsize=fontsize, **kw)
            texts.append(text)
    return texts


NUMPRIM2FONTSIZE = {i: 7 for i in range(4)}
NUMPRIM2FONTSIZE.update({i: 5 for i in range(4, 7)})
NUMPRIM2FONTSIZE.update({i: 3 for i in range(7, 10)})


def flexible_font_adjustment(length: int):
    if length >= 10:
        return 1
    else:
        return NUMPRIM2FONTSIZE[length]


def do_alignment(P1, P2):
    """test fwd and rev alignment"""
    correspondences = [
        ("R", "A"),
        ("R", "G"),
        ("Y", "C"),
        ("Y", "T"),
        ("M", "A"),
        ("M", "C"),
        ("K", "G"),
        ("K", "T"),
        ("S", "G"),
        ("S", "C"),
        ("W", "A"),
        ("W", "T"),
        ("B", "C"),
        ("B", "G"),
        ("B", "T"),
        ("D", "A"),
        ("D", "G"),
        ("D", "T"),
        ("H", "A"),
        ("H", "C"),
        ("H", "T"),
        ("V", "A"),
        ("V", "C"),
        ("V", "G"),
        ("N", "A"),
        ("N", "C"),
        ("N", "G"),
        ("N", "T"),
    ]  # Set correspondences between degenerated bases for edlib to know what to retain as equivalent while aligning: option is suggested by the PyPi page of the edlib package
    aln1 = edlib.align(
        P1, P2, task="path", mode="NW", additionalEqualities=correspondences
    )
    nice = edlib.getNiceAlignment(
        aln1, P1, P2
    )
    # print("\n".join(nice.values()))
    return aln1["editDistance"], "/".join(nice.values())


def read_primer_table(csv):
    """Read primers table, provided as a csv file whose separatore MUST be comma and whose fields MUST be named and ordered as follows: ID,FWD,REV
    ID: contains the demultiplexing unit to which the primers are referred
    FWD and REV: idicate respectively the forward and reverse primer sequences"""
    print("Reading primers table...", file=sys.stderr)
    if os.path.splitext(csv)[1] == ".csv":
        try:
            primers = pd.read_csv(csv)  # read primers table thanks to pandas
            fp = primers["FWD"]  # Forward sequences
            rp = primers["REV"]  # Reverse sequences
            ids = primers["ID"]  # Demultiplexing unit ID
            primer_dict = d(list)
            for i in range(len(ids)):
                # Assign to each demultiplexing unit ID its reverse and forward primers
                primer_dict[ids[i]] = [fp[i], rp[i]]
            print("Done", file=sys.stderr)
            return primer_dict
        except KeyError or ValueError:
            raise KeyError(
                "Fields of the csv do not comply with the requirements")
    else:
        raise ValueError("Input file is the wrong format")


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def calculate_pairwise_editdist(primersdict: dict, outputfile: str):
    primers = {f"{idp}_fwd": primersdict[idp][0]
               for idp in list(primersdict.keys())}
    primers.update(
        {f"{idp}_rev": primersdict[idp][1] for idp in list(primersdict.keys())}
    )
    primers.update(
        {f"rc_{idp}": reverse_complement(
            primers[idp]) for idp in list(primers.keys())}
    )
    lens = [len(primers[idp]) for idp in list(primers.keys())]
    matrix_size = len(primers)
    matrix = np.zeros((matrix_size, matrix_size))
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = np.nan
    compared_couples = []
    df = {"P1": [], "P2": [], "EDITDIST": [], "P1/ALN/P2": []}
    trees = list(primers.keys())

    for i in range(len(trees)):
        for j in range(len(trees)):
            if i != j:
                if [trees[i], trees[j]] not in compared_couples:
                    compared_couples.append([trees[i], trees[j]])
                    compared_couples.append([trees[j], trees[i]])
                    dist, ALN = do_alignment(
                        primers[trees[i]], primers[trees[j]])
                    df["P1"].append(trees[i])
                    df["P2"].append(trees[j])
                    matrix[i, j] = dist / max(lens)
                    df["EDITDIST"].append(dist / max(lens))
                    df["P1/ALN/P2"].append(ALN)
                else:
                    continue
            else:
                continue
    df = pd.DataFrame.from_dict(df)
    df.to_csv(outputfile, index=False)
    if not manyp or (manyp and matrix_size / 4 < 15):
        fig, ax = plt.subplots()
        im, cbar = heatmap(
            matrix,
            trees,
            trees,
            ax=ax,
            cmap="OrRd",
            cbarlabel="Primers edit dist",
        )
        if not humanread or (humanread and matrix_size / 4 < 10):
            texts = annotate_heatmap(
                im, valfmt="{x:.2f}", fontsize=flexible_font_adjustment(matrix_size / 4)
            )
        if humanread and matrix_size / 4 > 10:
            warnings.warn(
                "--human_readable is enabled and you have more than 10 pairs of primers, heatmap won't be annotated",
                HumanReadableWarning,
            )
        fig.tight_layout()
        fig.savefig(outputfile.replace(".csv", ".png"), dpi=300)
        return round(min(df["EDITDIST"]), 4)
    else:
        warnings.warn(
            "--many_primers is enabled and you have more than 15 pairs of primers, heatmap won't be generated",
            ManyPrimersWarning,
        )


if __name__ == "__main__":
    primersdict = read_primer_table(prim)
    print(calculate_pairwise_editdist(primersdict, outf))

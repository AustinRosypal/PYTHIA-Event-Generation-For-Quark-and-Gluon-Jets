import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Bringing in the output, merged root file from the job
root_file = "/home/ajrosy/HEP/FinalProject736/pythiaOutput.root"

plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['figure.titlesize'] = 18

def compute_stats(values, edges):
    centers = 0.5 * (edges[:-1] + edges[1:])
    entries = np.sum(values)

    if entries == 0:
        return 0, 0, 0

    mean = np.sum(values * centers) / entries
    variance = np.sum(values * (centers - mean)**2) / entries
    std = np.sqrt(variance)

    return entries, mean, std

def plot_th1(
    ax,
    hist,
    title,
    xlabel="Value",
    ylabel="Counts",
    x_min=None
):
    values, edges = hist.to_numpy()
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    ax.bar(
        centers,
        values,
        width=widths,
        align="center",
        color="darkorange",
        edgecolor="grey"
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")

    if x_min is not None:
        ax.set_xlim(left=x_min)

    # Compute stats
    entries, mean, std = compute_stats(values, edges)

    # Stats box
    stats_text = (
        f"Entries = {entries:.0f}\n"
        f"Mean = {mean:.3f}\n"
        f"Std Dev = {std:.3f}"
    )

    ax.text(
        0.97, 0.95,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(facecolor="white", edgecolor="black")
    )

def plot_th2(
    ax,
    hist,
    title,
    xlabel="X",
    ylabel="Y",
    draw_rect=False
):
    values, xedges, yedges = hist.to_numpy()

    mesh = ax.pcolormesh(xedges, yedges, values.T)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.colorbar(mesh, ax=ax)

    if draw_rect:
        y_min = yedges[0]
        y_max = yedges[-1]

        rect = patches.Rectangle(
            (-3.5, y_min),
            2.5,
            y_max - y_min,
            linewidth=4,
            edgecolor="red",
            facecolor="none"
        )

        ax.add_patch(rect)


# -------------------------------
# Open ROOT file and plot
# -------------------------------
with uproot.open(root_file) as f:

    h_quarkMult = f["h_quarkMult"]
    h_gluonMult = f["h_gluonMult"]

    h_quarkAngDist = f["h_quarkGirth"]
    h_gluonAngDist = f["h_gluonGirth"]
    #h_quarkAngDist = f["h_quarkAngDist"]
    #h_gluonAngDist = f["h_gluonAngDist"]

    h_jetpt = f["h_jetpt"]

    h_quarkPhaseSpace = f["h_quarkPhaseSpace"]
    h_gluonPhaseSpace = f["h_gluonPhaseSpace"]

    # -------------------------------
    # Multiplicity
    # -------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    plot_th1(
        axes[0],
        h_quarkMult,
        "Quark Jet Charged Particle Multiplicity",
        xlabel="Charged Particle Multiplicity"
    )

    plot_th1(
        axes[1],
        h_gluonMult,
        "Gluon Jet Charged Particle Multiplicity",
        xlabel="Charged Particle Multiplicity"
    )

    fig.suptitle("Quark and Gluon Jet Multiplicity Comparison")
    plt.tight_layout()
    plt.show()

    # -------------------------------
    # Angular distributions / Girth
    # -------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    plot_th1(
        axes[0],
        h_quarkAngDist,
        "Quark Jet Angular Girth",
        xlabel="Girth"
    )

    plot_th1(
        axes[1],
        h_gluonAngDist,
        "Gluon Jet Angular Girth",
        xlabel="Girth"
    )

    fig.suptitle("Quark and Gluon Jet Angular Girth")
    plt.tight_layout()
    plt.show()

    # -------------------------------
    # Jet p_T
    # -------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    plot_th1(
        ax,
        h_jetpt,
        "Jet $p_T$",
        xlabel=r"$p_{T}$ [GeV/c]",
        x_min=5
    )

    plt.tight_layout()
    plt.show()

    # -------------------------------
    # Geometrical Phase-space histograms
    # -------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    plot_th2(
        axes[0],
        h_quarkPhaseSpace,
        "Quark Jet Geometrical Phase Space",
        xlabel=(r"$\eta$"),
        ylabel=(r"$\phi$"),
        draw_rect=True
    )

    plot_th2(
        axes[1],
        h_gluonPhaseSpace,
        "Gluon Jet Geometrical Phase Space",
        xlabel=(r"$\eta$"),
        ylabel=(r"$\phi$"),
        draw_rect=True
    )

    fig.suptitle("Quark and Gluon Jet Geometrical Phase Space Distributions")
    plt.tight_layout()
    plt.show()



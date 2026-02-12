import matplotlib.pyplot as plt
import seaborn as sns

def set_plot_style():
    plt.rcParams['figure.figsize'] = (5, 4)
    plt.rcParams['figure.dpi'] = 100

    sns.set_style("ticks")
    sns.set_theme(
        style="white",
        context="notebook",
        font_scale=1.2,
        rc={
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "font.size": 14,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 12,
            "mathtext.fontset": "cm",
            "xtick.bottom": True,
            "ytick.left": True,
            "axes.titleweight": "bold",
            "legend.frameon": False,
            "legend.framealpha": 1,
            "legend.facecolor": "white",
            "legend.edgecolor": "black",
            "lines.linewidth": 1.75,
            "axes.linewidth": 1.2,
            "grid.linewidth": 0.6,
            "grid.alpha": 0.8,
            "grid.linestyle": "--",
            "lines.markersize": 5,
        }
    )

    sns.set_palette("bright")

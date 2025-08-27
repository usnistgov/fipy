from matplotlib import pyplot as plt
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

color_map = {
    'no-pysparse': ('red', "s"),
    'trilinos': ('red', "s"),
    'petsc': ('blue', "o"),
    'petsc 3.18.5': ('blue', "o"),
    'petsc 3.20.2': ('turquoise', "H"),
    'scipy': ('green', "v"),
    'pysparse': ('orange', "^"),
    'pyamgx': ('cyan', "*"),
    'petsc-RCV': ('pink', "o")
}

# line_style = {
#     'a5f233aa7': ("", "full"),
#     '371d28468': ("", "none"),
#     "ef3db887e": ("", "full")
# }
line_style = {
    'a5f233aa7': ("-", "full"),
    '371d28468': ("--", "none"),
    "ef3db887e": ("-.", "full")
}

def plot_scaling(df, color_by_suite=True,
                 by=["package.solver", "solver_class", "preconditioner"],
                 xdata="tasks",
                 ydata="elapsed / s", 
                 ax=None,
                 linewidth=1., alpha=1., show_marker=True,
                 style="sigma"):

    if ax is None:
        fig, ax = plt.subplots(figsize=(8,6),
                               gridspec_kw={"right": 0.8})
    groups = df.groupby(by + [xdata])
    groups = groups.agg(data_count=(ydata, "count"),
                        data_mean=(ydata, "mean"),
                        data_min=(ydata, "min"),
                        data_max=(ydata, "max"),
                        data_std=(ydata, "std")).reset_index()
    groups = groups.groupby(by)

    for keys, group in groups:
        suite = keys[0]
        fipy_rev = keys[1]
        if color_by_suite:
            color, marker = color_map[suite]
            if not show_marker:
                marker = ""
            line, fill = line_style[fipy_rev]
        group.plot(xdata, "data_mean",
                   ax=ax, label=None, color=color, marker=marker,
                   linestyle=line, fillstyle=fill,
                   linewidth=linewidth, alpha=alpha)

        if style == "minmax":
            # plot range
            if len(group) == 1:
                ax.errorbar(group[xdata], group["data_mean"],
                            yerr=[group["data_mean"] - group["data_min"],
                                  group["data_max"] - group["data_mean"]], color=color)
            else:
                ax.fill_between(group[xdata],
                                group["data_min"],
                                group["data_max"],
                                color=color,
                                alpha=0.1)
        elif style == "sigma":
            # plot uncertainty
            err = group["data_std"] / np.sqrt(group["data_count"])
            if len(group) == 1:
                ax.errorbar(group[xdata], group["data_mean"], yerr=err, color=color)
            else:
                ax.fill_between(group[xdata],
                                group["data_mean"] - err,
                                group["data_mean"] + err,
                                color=color,
                                alpha=0.1)

    return ax

def suite_legend(suites, ax, color_by_suite=True, loc="upper left", linestyle="-"):
    if color_by_suite:
        legend_elements = [Line2D([0], [0],
                                  color=c, label=s, marker=m,
                                  linestyle=linestyle)
                           for s, (c, m) in color_map.items()
                           if s in suites.unique()]
        ax.legend(handles=legend_elements, loc=loc)
    else:
        # only label converged lines
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::2], labels[::2],  loc=loc)

def revision_legend(ax, loc="lower right"):
    legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle="--", label="FiPy 3.4.4 (371d28468)"),
                       Line2D([0], [0], color="black", marker="o", label="FiPy 4.0 (a5f233aa7)"),
                       Patch(facecolor="black", edgecolor=None, alpha=0.1, label="$\pm 1$ standard deviation"),
                       Line2D([0], [0], color="black", linestyle="-", linewidth=1, alpha=0.1, label=r"$\sim N$"),
                       Line2D([0], [0], color="black", linestyle=":", linewidth=1, alpha=0.1, label=r"$N \log N$")]
    leg = Legend(ax, handles=legend_elements, labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)", "$\pm 1$ standard deviation", r"$\sim N$", r"$\sim N\, \log N$"],
                 loc=loc) #, frameon=False)
    ax.add_artist(leg)

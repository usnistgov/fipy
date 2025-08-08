from __future__ import unicode_literals
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.patches import Patch

import pandas as pd
from scipy.io import mmread

def plot_all(df, output, color_by_suite=True,
             by=["package.solver", "solver_class", "preconditioner"],
             xdata="tasks", xlabel="tasks",
             ydata="elapsed / s", ylabel="elapsed time", title=None,
             xmin=None, xmax=None, ymin=None, ymax=None, ax=None,
             linewidth=1., alpha=1., show_marker=True, legends=True,
             style="sigma", logx=True, logy=True):
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
    line_style = {
        'a5f233aa7': ("-", "full"),
        '371d28468': ("--", "none"),
        "ef3db887e": ("-.", "full")
    }

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
            group.plot(xdata, "data_mean", logx=logx, logy=logy,
                       ax=ax, label=None, color=color, marker=marker, fillstyle=fill,
                       linestyle=line, linewidth=linewidth, alpha=alpha)
        else:
            group.plot(xdata, "data_mean", logx=logx, logy=logy,
                       ax=ax, label=None, color=color, marker=marker, fillstyle=fill,
                       linestyle=line, linewidth=linewidth, alpha=alpha)

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

    if logx & logy:
        N = np.logspace(3, 6, 100)
        ax.loglog(N, N * 1e-3, linewidth=0.5, color="black")
        ax.text(3.3e3, 4.5e0, r"$\sim N$", rotation=45)
        ax.loglog(N, N * np.log(N) * 1e-8, linewidth=0.5, color="black")
        ax.text(3.3e4, 2.3e-3, r"$\sim N\, \ln N$", rotation=45)

    if "solver" in legends:
        if color_by_suite:
            legend_elements = [Line2D([0], [0], color=c, label=s, marker=m)
                               for s, (c, m) in color_map.items()
                               if s in df["package.solver"].unique()]
            ax.legend(handles=legend_elements, loc="upper left")
        else:
            # only label converged lines
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::2], labels[::2],  loc="upper left")
    else:
        ax.get_legend().remove()

    if "revs" in legends:
        legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle="--"),
                           Line2D([0], [0], color="black", marker="o"),
                           Patch(facecolor="black", edgecolor=None, alpha=0.1)]
        leg = Legend(ax, handles=legend_elements, labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)", "minimum to maximum"],
                     loc='lower left') #, frameon=False)
        ax.add_artist(leg)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    if title is not None:
        ax.set_title(title)

    if output is not None:
        plt.savefig(output)

    return ax

if __name__ == "__main__":
    all = pd.read_json("serial_scaling.json")
    all = all[all["fipy_rev"].isin(["371d28468", "a5f233aa7"])]

    all.loc[all["package.solver"] == "petsc",
            "package.solver"] = all.apply(lambda r: f"petsc {r['package.petsc4py']}", axis=1)
    all.loc[all["package.solver"] == "no-pysparse",
            "package.solver"] = "trilinos"

    all = all.copy() # defrag

    all["prepare2solve"] = all["prepare_seconds"] / all["solve_seconds"]
    all["prepare2elapsed"] = all["prepare_seconds"] / all["elapsed_seconds"]

    fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(6, 9))

    for suite, ax, subfig in [["pysparse", axs[0,0], "a"],
                              ["scipy", axs[0,1], "b"],
                              ["petsc 3.18.5", axs[1,0], "c"],
                              ["petsc 3.20.2", axs[1,1], "d"],
                              ["trilinos", axs[2,0], "e"]]:
        suite_data = all[all["package.solver"] == suite]
        plot_all(suite_data, None,
                 by=["package.solver", "fipy_rev", "solver_class", "preconditioner"],
                 xdata="numberOfElements", xlabel="number of cells",
                 ydata="prepare2elapsed", ylabel="prepare time / elapsed time",
                 ymin=0, ymax=1, ax=ax, title=f"({subfig}) {suite}",
                 linewidth=1, alpha=0.1, style="none", show_marker=False,
                 logy=False, legends="")
        plot_all(suite_data[((suite_data["solver_class"] == "LinearGMRESSolver")
                             & (suite_data["preconditioner"] == "JacobiPreconditioner"))], None,
                 by=["package.solver", "fipy_rev"],
                 xdata="numberOfElements", xlabel="number of cells",
                 ydata="prepare2elapsed", ylabel="prepare time / elapsed time",
                 ymin=0, ymax=1, ax=ax, title=f"({subfig}) {suite}",
                 linewidth=2, style="none", show_marker=False,
                 logy=False, legends="")

    plt.tight_layout()
    plt.show()

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

    if legends:
        if color_by_suite:
            legend_elements = [Line2D([0], [0], color=c, label=s, marker=m)
                               for s, (c, m) in color_map.items()
                               if s in df["package.solver"].unique()]
            ax.legend(handles=legend_elements, loc="upper left")
        else:
            # only label converged lines
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::2], labels[::2],  loc="upper left")

        legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle="--", label="FiPy 3.4.4 (371d28468)"),
                           Line2D([0], [0], color="black", marker="o", label="FiPy 4.0 (a5f233aa7)"),
                           Patch(facecolor="black", edgecolor=None, alpha=0.1, label="$\pm 1$ standard deviation"),
                           Line2D([0], [0], color="black", linestyle="-", linewidth=1, alpha=0.1, label=r"$\sim N$"),
                           Line2D([0], [0], color="black", linestyle=":", linewidth=1, alpha=0.1, label=r"$N \log N$")]
        leg = Legend(ax, handles=legend_elements, labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)", "$\pm 1$ standard deviation", r"$\sim N$", r"$\sim N\, \log N$"],
                     loc='lower right') #, frameon=False)
        ax.add_artist(leg)
    else:
        ax.get_legend().remove()

    ax.autoscale(False)
    if logx & logy:
        N = np.logspace(0, 7, 100)
        for n in range(-10, 3, 2):
            ax.axline((1, 10**n), (10, 10**(n+1)), color="black", linewidth=1, alpha=0.1)
        for n in range(-14, 4, 2):
            ax.loglog(N, N * np.log(N) * 10**n, linewidth=1, color="black", alpha=0.1, linestyle=":")

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

    fig = plt.figure(figsize=(10, 10))
    outer_grid = fig.add_gridspec(2, 2)
    axs = outer_grid.subplots()
    
    axs[0,0].sharey(axs[0,1])
    axs[0,0].sharex(axs[1,0])

    plot_all(all, None, by=["package.solver", "fipy_rev"],
             xdata="numberOfElements", xlabel="number of cells",
             ydata="elapsed_seconds", ylabel="time / s",
             ymin=1e-4, ymax=1e2, ax=axs[0, 0], title="(a) elapsed time")

    plot_all(all, None, by=["package.solver", "fipy_rev"],
             xdata="numberOfElements", xlabel="number of cells",
             ydata="prepare_seconds", ylabel="time / s",
             ymin=1e-4, ymax=1e2, ax=axs[0,1], title="(b) prepare time", legends=False)

    plot_all(all, None, by=["package.solver", "fipy_rev"],
             xdata="numberOfElements", xlabel="number of cells",
             ydata="solve_seconds", ylabel="time / s",
             ymin=1e-4, ymax=1e2, ax=axs[1,0], title="(c) solve time", legends=False)

#     plot_all(all, None,
#              by=["package.solver", "fipy_rev", "solver_class", "preconditioner"],
#              xdata="numberOfElements", xlabel="number of cells",
#              ydata="solve_seconds", ylabel="time / s",
#              ymin=1e-4, ymax=1e2, linewidth=0.2, alpha=0.2, style="none", show_marker=False,
#              ax=axs[1,0], title="(c) solve time", legends=False)
#     plot_all(all[((all["solver_class"] == "LinearGMRESSolver")
#                   & (all["preconditioner"] == "JacobiPreconditioner"))], None,
#              by=["package.solver", "fipy_rev"],
#              xdata="numberOfElements", xlabel="number of cells",
#              ydata="solve_seconds", ylabel="time / s",
#              ymin=1e-4, ymax=1e2, style="none", linewidth=2,
#              ax=axs[1,0], title="(c) solve time", legends=False)

#     legend_elements = [Line2D([0], [0], color="black", marker="", linewidth=0.2, alpha=0.2),
#                        Line2D([0], [0], color="black", marker="o")]
#     leg = Legend(axs[1,0], handles=legend_elements, labels=["each solver & preconditioner", "LinearGMRESSolver & JacobiPreconditioner"],
#                  loc='upper left') #, frameon=False)
#     axs[1,0].add_artist(leg)
             
    inner_fig = fig.add_subfigure(outer_grid[1, 1])
    inner_grid = inner_fig.add_gridspec(10, 10)
    
    sparse_fig = inner_fig.add_subfigure(inner_grid[1:5, :6])
    sparse_fig.suptitle("sparse linear system")
    
    sparse_grid = sparse_fig.add_gridspec(16, 24)
    
    L_ax = sparse_fig.add_subplot(sparse_grid[:, :16])
    c_ax = sparse_fig.add_subplot(sparse_grid[:, -5])
    b_ax = sparse_fig.add_subplot(sparse_grid[:, 16:-5])
    soln_ax = inner_fig.add_subplot(inner_grid[-5:-1, -5:-1])
    
    coo = mmread("serial_scaling.mtx")
    rhs = np.load("serial_scaling.rhs.npz")

    L = coo.todense()
    b = rhs["rhs"][:, np.newaxis]
    loc = None
    
    # t = SymmetricalLogTransform(base=10, linthresh=1, linscale=1)
    # L = t.transform(coo.todense().flat).reshape((36, 36))
    # b = t.transform(rhs["rhs"]).reshape((36, 1))
    # loc = SymmetricalLogLocator(base=10, linthresh=1)

    zRange = max(abs(L).max(), abs(b).max())

    L_ax.matshow(L, cmap="RdBu", vmin=-zRange, vmax=zRange)
    b_ax.imshow(b, cmap="RdBu", vmin=-zRange, vmax=zRange, aspect=0.5)
    
    norm = Normalize(vmin=-zRange, vmax=zRange)
    ColorbarBase(ax=c_ax, cmap="RdBu", norm=norm, orientation='vertical',
                 format=None, ticks=loc)
                 
    L_ax.text(0.5, -0.1, "L",
              transform=L_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')
    b_ax.text(0.5, -0.1, "b",
              transform=b_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')

    for ax in [L_ax, b_ax, soln_ax]:
        ax.set(xticks=[], yticks=[])
    
    axs[1, 1].axis('off')

    data = np.load("serial_initial.npz")
    soln_ax.imshow(0.7 - 0.4 * data["phi"])
    soln_ax.set_title("initial condition")

    plt.tight_layout()
    plt.show()

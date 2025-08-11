from __future__ import unicode_literals

from functools import partial
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.patches import Patch
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.ticker import SymmetricalLogLocator, LogFormatter, LogFormatterMathtext
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import pandas as pd
from scipy.io import mmread
from scipy.optimize import curve_fit


def plot_all(df, output, color_by_suite=True,
             by=["package.solver", "solver_class", "preconditioner"],
             xdata="tasks", xlabel="tasks",
             ydata="elapsed / s", ylabel="elapsed time", title=None,
             xmin=None, xmax=None, ymin=None, ymax=None, ax=None,
             linewidth=1., show_marker=True,
             style="sigma"):
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
        'a5f233aa7': ("", "full"),
        '371d28468': ("", "none"),
        "ef3db887e": ("", "full")
    }

    # matplotlib.use('Agg')

    # plt.figure()
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
            group.plot(xdata, "data_mean", loglog=True,
                       ax=ax, label=None, color=color, marker=marker, linestyle=line, fillstyle=fill,
                       linewidth=linewidth)
        else:
            group.plot(xdata, "data_mean", loglog=True,
                       ax=ax, label=None, color=color, marker=marker, linestyle=line, fillstyle=fill,
                       linewidth=linewidth)

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

    if color_by_suite:
        legend_elements = [Line2D([0], [0], color=c, label=s, marker=m, linestyle="")
                           for s, (c, m) in color_map.items()
                           if s in df["suite"].unique()]
        ax.legend(handles=legend_elements, loc="upper left")
    else:
        # only label converged lines
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::2], labels[::2],  loc="upper left")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    if title is not None:
        ax.set_title(title)

    # plt.show()

    if output is not None:
        plt.savefig(output)

    return ax

def amdahl(p, sigma, baseline=1):
    return baseline * p / (1 + sigma * (p - 1))

def USL(p, sigma, kappa, baseline=1):
    return baseline * p / (1 + sigma * (p - 1) + kappa * p * (p - 1))

if __name__ == "__main__":
    all = pd.read_json("parallel_scaling.json")

    all.loc[all["suite"] == "no-pysparse", "suite"] = "trilinos"

    by = ["suite", "fipy_rev", "nx"]
    data_set="elapsed / s"
    groups = all.groupby(by + ["tasks"])
    groups = groups.agg(data_mean=(data_set, "mean")).reset_index()
    groups = groups.groupby(by)
    for nx in all["nx"].unique():
        baseline = groups.get_group(("pysparse", "371d28468", nx)).iloc[0]["data_mean"]
        all.loc[all["nx"] == nx, "speedup"] = baseline / all[data_set]
        
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    plot_all(all[all["nx"] == 2048], None, by=["suite", "fipy_rev"],
             xdata="tasks", xlabel="tasks", ax=ax,
             ydata="speedup", ylabel="speedup / ($t_{PySparse} / t_N$)",
             ymin=1e-1, ymax=3e1)

#     d.groupby("tasks")["speedup"].mean().iloc[0]
#     (d.groupby("tasks")["speedup"].mean().iloc[0] / 2)

    fit_tasks = np.logspace(0, 2, 100)

    groups = all[(all["nx"] == 2048)].groupby(by=["suite", "fipy_rev"])
    
    for suite, rev, n, color in [["petsc", "a5f233aa7", 1, "blue"],
                                 ["trilinos", "a5f233aa7", 0, "red"],
                                 ["petsc", "371d28468", 0, "blue"],
                                 ["trilinos", "371d28468", 0, "red"]]:
        print(f"{suite} {rev}")
        
        d = groups.get_group((suite, rev))

        d = d[d["tasks"] >= 2**n]
        
        baseline = d.groupby("tasks")["speedup"].mean().iloc[0] / 2**n

#         ax.plot(d["tasks"], d["speedup"], marker="x", linestyle="", color=color)

        amdahl_bl = partial(amdahl, baseline=baseline)
        popt, pcov = curve_fit(amdahl_bl, d["tasks"], d["speedup"])
        ax.plot(fit_tasks, amdahl_bl(fit_tasks, *popt), color=color, linestyle="--")
        print(f"Amdahl")
        print(f"sigma = {popt[0]} +/- {np.sqrt(pcov[0][0])}")

        USL_bl = partial(USL, baseline=baseline)
        popt, pcov = curve_fit(USL_bl, d["tasks"], d["speedup"])
        ax.plot(fit_tasks, USL_bl(fit_tasks, *popt), color=color, linestyle="-")
        print(f"Gunther")
        print(f"sigma={popt[0]} +/- {np.sqrt(pcov[0][0])}")
        print(f"kappa={popt[1]} +/- {np.sqrt(pcov[1][1])}")

    legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle=""),
                       Line2D([0], [0], color="black", marker="o", linestyle=""),
                       Patch(facecolor="black", edgecolor=None, alpha=0.1),
                       Line2D([0], [0], color="black", marker="", linestyle="--"),
                       Line2D([0], [0], color="black", marker="", linestyle="-"),
                       Line2D([0], [0], color="black", linestyle="-", linewidth=1, alpha=0.1)]
    leg = Legend(ax, handles=legend_elements,
                 labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)",
                         "$\pm 1$ standard deviation",
                         "Amdahl's Law", "Gunther's Law", r"$\sim N$"],
                 loc=(0.7, 0.35)) #, frameon=False)
    ax.add_artist(leg);
                 
    for n in range(-3, 2):
        ax.axline((1, 10**n), (10, 10**(n+1)), color="black", linewidth=1, alpha=0.1)

    outer_grid = fig.add_gridspec(3, 3)
    inner_fig = fig.add_subfigure(outer_grid[2, 1])
    inner_grid = inner_fig.add_gridspec(10, 10)
    
    sparse_fig = inner_fig.add_subfigure(inner_grid[:6, :9])
    sparse_fig.suptitle("sparse linear system")
    
    sparse_grid = sparse_fig.add_gridspec(16, 24)
    
    L_ax = sparse_fig.add_subplot(sparse_grid[:, :16])
    c_ax = sparse_fig.add_subplot(sparse_grid[:, -5])
    b_ax = sparse_fig.add_subplot(sparse_grid[:, 16:-5])
    
    coo = mmread("parallel_scaling.mtx")
    rhs = np.load("parallel_scaling.rhs.npz")

    L = coo.todense()
    b = rhs["rhs"][:, np.newaxis]
    loc = None
    fmt = None
    
    t = SymmetricalLogTransform(base=10, linthresh=0.0001, linscale=10)
    L = t.transform(L.flat).reshape((32, 32))
    b = t.transform(b).reshape((32, 1))
    loc = SymmetricalLogLocator(transform=t)
    fmt = LogFormatterMathtext(linthresh=0.0001)

    zRange = max(abs(L).max(), abs(b).max())

    L_ax.matshow(L, cmap="RdBu", vmin=-zRange, vmax=zRange)
    b_ax.imshow(b, cmap="RdBu", vmin=-zRange, vmax=zRange, aspect=0.5)
    
    norm = Normalize(vmin=-zRange, vmax=zRange)
    ColorbarBase(ax=c_ax, cmap="RdBu", norm=norm, orientation='vertical',
                 format=fmt, ticks=loc)
                 
    L_ax.text(0.5, -0.1, "L",
              transform=L_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')
    b_ax.text(0.5, -0.1, "b",
              transform=b_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')

              
    soln_ax = fig.add_subplot(outer_grid[2, 2], position=[0.645, 0.11, 0.23, 0.23])

    data = np.load("parallel_initial.npz")
    soln_ax.imshow(data["eta"])
    soln_ax.set_title("initial condition")
    
    for ax in [L_ax, b_ax, soln_ax]:
        ax.set(xticks=[], yticks=[])
    
#     axs[1, 1].axis('off')
            
#     popt, pcov = curve_fit(USL, d["tasks"], d["speedup"], (d.groupby("tasks")["speedup"].mean().iloc[0] / 2,))
#     ax.plot(fit_tasks, USL(fit_tasks, *popt), color="blue", linestyle="-")
    
#     popt, pcov = curve_fit(amdahl, d["tasks"], d["speedup"])
#     fit_tasks = np.logspace(0, 2, 100)
#     ax.plot(fit_tasks, amdahl(fit_tasks, *popt), color="blue", linestyle="--")

    
# #     legend_elements = [Line2D([0], [0], color="black", marker="", linewidth=0.2, alpha=0.2),
# #                        Line2D([0], [0], color="black", marker="o")]
# #     leg = Legend(axs[1,0], handles=legend_elements, labels=["each solver & preconditioner", "LinearGMRESSolver & JacobiPreconditioner"],
# #                  loc='upper left') #, frameon=False)
# #     axs[1,0].add_artist(leg)
#              
#     inner_fig = fig.add_subfigure(outer_grid[1, 1])
#     inner_grid = inner_fig.add_gridspec(10, 10)
#     
#     sparse_fig = inner_fig.add_subfigure(inner_grid[1:5, :6])
#     sparse_fig.suptitle("sparse linear system")
#     
#     sparse_grid = sparse_fig.add_gridspec(16, 24)
#     
#     L_ax = sparse_fig.add_subplot(sparse_grid[:, :16])
#     c_ax = sparse_fig.add_subplot(sparse_grid[:, 16])
#     b_ax = sparse_fig.add_subplot(sparse_grid[:, 17:])
#     soln_ax = inner_fig.add_subplot(inner_grid[-5:-1, -5:-1])
#     
#     coo = mmread("serial_scaling.mtx")
#     rhs = np.load("serial_scaling.rhs.npz")
# 
#     L = coo.todense()
#     b = rhs["rhs"][:, np.newaxis]
#     loc = None
#     
#     # t = SymmetricalLogTransform(base=10, linthresh=1, linscale=1)
#     # L = t.transform(coo.todense().flat).reshape((36, 36))
#     # b = t.transform(rhs["rhs"]).reshape((36, 1))
#     # loc = SymmetricalLogLocator(base=10, linthresh=1)
# 
#     zRange = max(abs(L).max(), abs(b).max())
# 
#     L_ax.matshow(L, cmap="RdBu", vmin=-zRange, vmax=zRange)
#     b_ax.imshow(b, cmap="RdBu", vmin=-zRange, vmax=zRange, aspect=0.5)
#     
#     norm = Normalize(vmin=-zRange, vmax=zRange)
#     ColorbarBase(ax=c_ax, cmap="RdBu", norm=norm, orientation='vertical',
#                  format=None, ticks=loc)
#                  
#     L_ax.text(0.5, -0.1, "L",
#               transform=L_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')
#     b_ax.text(0.5, -0.1, "b",
#               transform=b_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')
# 
#     for ax in [L_ax, b_ax, soln_ax]:
#         ax.set(xticks=[], yticks=[])
#     
#     axs[1, 1].axis('off')
# 
#     data = np.load("serial_initial.npz")
#     soln_ax.imshow(0.7 - 0.4 * data["phi"])
#     soln_ax.set_title("initial condition")

#     plt.tight_layout()
    plt.show()

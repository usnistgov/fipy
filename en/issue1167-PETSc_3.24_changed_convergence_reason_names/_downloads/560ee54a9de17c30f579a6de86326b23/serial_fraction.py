import numpy as np

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend import Legend

import pandas as pd

from plot_scaling import ScalePlot

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

#     fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(9, 6))
    fig = plt.figure(figsize=(9, 6))
    grid = fig.add_gridspec(2, 3)
    axs = grid.subplots()

    axs[0,0].sharey(axs[0,1])
    axs[0,1].sharey(axs[0,2])
    axs[1,0].sharey(axs[1,1])
    axs[0,0].sharex(axs[1,0])
    axs[0,1].sharex(axs[1,1])

    for suite, ax, letter in [["pysparse", axs[0,0], "a"],
                              ["scipy", axs[1,0], "b"],
                              ["petsc 3.18.5", axs[0,1], "c"],
                              ["petsc 3.20.2", axs[1,1], "d"],
                              ["trilinos", axs[0,2], "e"]]:
        suite_data = all[all["package.solver"] == suite]

        scaling = ScalePlot(ax=ax)
        scaling.plot(suite_data,
                     by=["package.solver", "fipy_rev", "solver_class", "preconditioner"],
                     xdata="numberOfElements",
                     ydata="prepare2elapsed",
                     linewidth=1, alpha=0.1, show_marker=False)

        scaling.plot(suite_data[(suite_data["solver_class"].isin(["LinearCGSolver", "LinearPCGSolver"])
                                 & (suite_data["preconditioner"] == "unpreconditioned"))],
                     by=["package.solver", "fipy_rev"],
                     xdata="numberOfElements",
                     ydata="prepare2elapsed",
                     linewidth=2, show_marker=False)

        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel("number of cells")
        ax.set_ylabel("prepare time / elapsed time")
        ax.set_xlim(xmin=None, xmax=None)
        ax.set_ylim(ymin=0, ymax=1)

        ax.set_title(f"({letter}) {suite}")
        ax.get_legend().remove()

    axs[1, 2].axis('off')

    legend_elements = [Line2D([0], [0], color="black", marker="", linestyle="--"),
                       Line2D([0], [0], color="black", marker="")]
    leg = Legend(axs[1, 2],
                 handles=legend_elements,
                 labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)"],
                 loc='lower left') #, frameon=False)
    axs[1, 2].add_artist(leg)

    legend_elements = [Line2D([0], [0], color="black", marker="", linewidth=1, alpha=0.1),
                       Line2D([0], [0], color="black", marker="", linewidth=2)]
    leg = Legend(axs[1, 2],
                 handles=legend_elements,
                 labels=["each solver & precoditioner",
                         "LinearCGSolver\nunpreconditioned"],
                 loc='upper left') #, frameon=False)
    axs[1, 2].add_artist(leg)

    plt.tight_layout()
    plt.show()

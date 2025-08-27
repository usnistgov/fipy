from __future__ import unicode_literals
import numpy as np

from matplotlib import pyplot as plt

import pandas as pd
from scipy.io import mmread

from plot_scaling import plot_scaling, suite_legend, revision_legend
from plot_sparse import plot_sparse

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

    plot_scaling(all, by=["package.solver", "fipy_rev"],
                 xdata="numberOfElements",
                 ydata="elapsed_seconds",
                 ax=axs[0, 0])

    axs[0,0].set_title("(a) elapsed time")

    plot_scaling(all, by=["package.solver", "fipy_rev"],
                 xdata="numberOfElements",
                 ydata="prepare_seconds",
                 ax=axs[0,1])

    axs[0,1].set_title("(b) prepare time")

    plot_scaling(all, by=["package.solver", "fipy_rev"],
                 xdata="numberOfElements",
                 ydata="solve_seconds",
                 ax=axs[1,0])

    axs[1,0].set_title("(c) solve time")

    for ax in [axs[0, 0], axs[0, 1], axs[1, 0]]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("number of cells")
        ax.set_ylabel("time / s")
        ax.set_xlim(xmin=None, xmax=None)
        ax.set_ylim(ymin=1e-4, ymax=1e2)

        N = np.logspace(0, 7, 100)
        for n in range(-10, 3, 2):
            ax.axline((1, 10**n), (10, 10**(n+1)), color="black", linewidth=1, alpha=0.1)
        for n in range(-14, 4, 2):
            ax.loglog(N, N * np.log(N) * 10**n, linewidth=1, color="black", alpha=0.1, linestyle=":")

    suite_legend(suites=all["package.solver"], ax=axs[0,0])
    revision_legend(ax=axs[0,0])

    for ax in [axs[0,1], axs[1,0]]:
        ax.get_legend().remove()

    inner_fig = fig.add_subfigure(outer_grid[1, 1])
    inner_grid = inner_fig.add_gridspec(10, 10)
    
    sparse_fig = inner_fig.add_subfigure(inner_grid[1:5, :6])
    coo = mmread("serial_scaling.mtx")
    rhs = np.load("serial_scaling.rhs.npz")["rhs"]
    plot_sparse(coo, rhs, fig=sparse_fig)
    
    axs[1, 1].axis('off')

    data = np.load("serial_initial.npz")
    soln_ax = inner_fig.add_subplot(inner_grid[-5:-1, -5:-1])
    soln_ax.imshow(0.7 - 0.4 * data["phi"])
    soln_ax.set_title("initial condition")
    soln_ax.set(xticks=[], yticks=[])

    plt.tight_layout()
    plt.show()

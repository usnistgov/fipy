from matplotlib import pyplot as plt
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

class ScalePlot:
    def __init__(self,
                 ax=None,
                 color_map={
                     'no-pysparse': ('red', "s"),
                     'trilinos': ('red', "s"),
                     'petsc': ('blue', "o"),
                     'petsc 3.18.5': ('blue', "o"),
                     'petsc 3.20.2': ('turquoise', "H"),
                     'scipy': ('green', "v"),
                     'pysparse': ('orange', "^"),
                     'pyamgx': ('cyan', "*"),
                     'petsc-RCV': ('pink', "o")
                 },
                 line_style = {
                     'a5f233aa7': ("-", "full"),
                     '371d28468': ("--", "none"),
                     "ef3db887e": ("-.", "full")
                 },
                 color_by_suite=True):
        self.color_map = color_map
        self.line_style = line_style
        self.color_by_suite = color_by_suite

        if ax is None:
            _, self.ax = plt.subplots(figsize=(8,6),
                                      gridspec_kw={"right": 0.8})
        else:
            self.ax = ax

    def _plot_dispersion(self, x, data_mean, data_min, data_max, data_err, color):
        if len(x) == 1:
            self.ax.errorbar(x, data_mean,
                             yerr=data_err,
                             color=color)
        else:
            self.ax.fill_between(x,
                                 data_min,
                                 data_max,
                                 color=color,
                                 alpha=0.1)

    def plot_dispersion(self, group, xdata, color):
        pass

    def plot(self, df,
             by=["package.solver", "solver_class", "preconditioner"],
             xdata="tasks",
             ydata="elapsed / s",
             linewidth=1., alpha=1., show_marker=True):

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
            if self.color_by_suite:
                color, marker = self.color_map[suite]
                if not show_marker:
                    marker = ""
                line, fill = self.line_style[fipy_rev]
            group.plot(xdata, "data_mean",
                       ax=self.ax, label=None, color=color, marker=marker,
                       linestyle=line, fillstyle=fill,
                       linewidth=linewidth, alpha=alpha)

            self.plot_dispersion(group, xdata, color)

    def suite_legend(self, suites, loc="upper left", linestyle="-"):
        if self.color_by_suite:
            legend_elements = [Line2D([0], [0],
                                      color=c, label=s, marker=m,
                                      linestyle=linestyle)
                               for s, (c, m) in self.color_map.items()
                               if s in suites.unique()]
            self.ax.legend(handles=legend_elements, loc=loc)
        else:
            # only label converged lines
            handles, labels = self.ax.get_legend_handles_labels()
            self.ax.legend(handles[::2], labels[::2],  loc=loc)

    def revision_legend(self, loc="lower right"):
        legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle="--", label="FiPy 3.4.4 (371d28468)"),
                           Line2D([0], [0], color="black", marker="o", label="FiPy 4.0 (a5f233aa7)"),
                           Patch(facecolor="black", edgecolor=None, alpha=0.1, label="$\pm 1$ standard deviation"),
                           Line2D([0], [0], color="black", linestyle="-", linewidth=1, alpha=0.1, label=r"$\sim N$"),
                           Line2D([0], [0], color="black", linestyle=":", linewidth=1, alpha=0.1, label=r"$N \log N$")]
        leg = Legend(self.ax, handles=legend_elements, labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)", "$\pm 1$ standard deviation", r"$\sim N$", r"$\sim N\, \log N$"],
                     loc=loc) #, frameon=False)
        self.ax.add_artist(leg)

class MinMaxScalePlot(ScalePlot):
    def plot_dispersion(self, group, xdata, color):
        # plot range
        self._plot_dispersion(x=group[xdata],
                              data_mean=group["data_mean"],
                              data_min=group["data_min"],
                              data_max=group["data_max"],
                              data_err=[group["data_mean"] - group["data_min"],
                                        group["data_max"] - group["data_mean"]],
                              color=color)

class UncertaintyScalePlot(ScalePlot):
    def plot_dispersion(self, group, xdata, color):
        # plot uncertainty
        err = group["data_std"] / np.sqrt(group["data_count"])
        self._plot_dispersion(x=group[xdata],
                              data_mean=group["data_mean"],
                              data_min=group["data_mean"] - err,
                              data_max=group["data_mean"] + err,
                              data_err=err,
                              color=color)

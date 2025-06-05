"""
1D Viewer that calculates and plots the position of the phase field
interface.
"""

from fipy import Matplotlib1DViewer, numerix
from scipy.optimize import leastsq
from matplotlib import pyplot as plt

class PhaseViewer(Matplotlib1DViewer):
    def __init__(self, phase, C, sharp, elapsed, L, deltaA, title=None,
                 tmin=None, tmax=None, **kwlimits):
        self.phase = phase
        self.x = self.phase.mesh.cellCenters[0]
        self.elapsed = elapsed
        self.L = L
        self.deltaA = deltaA

        fig = plt.figure(constrained_layout=True, figsize=[4, 6])
        gs = fig.add_gridspec(ncols=1, nrows=3)
        self.viewax = fig.add_subplot(gs[1:, 0])
        self.viewax.set_xlabel(r"$x / \mathrm{cm}$")
        self.posax = fig.add_subplot(gs[0, 0], sharex=self.viewax)
        self.posax.set_ylabel(r"$t / \mathrm{s}$")
        self.posax.label_outer()
        self.posax.set_ylim(tmin, tmax)

        self.times = []
        self.positions = []
        self.line, = self.posax.semilogy([0.], [0.],
                                         color="blue")
        self.bug, = self.posax.semilogy([], [],
                                        linestyle="", marker="o", color="red")

        super(PhaseViewer, self).__init__(vars=(phase, C, sharp),
                                          axes=self.viewax,
                                          **kwlimits)

        self.timer = self.viewax.text(0.00025, 0.2, "t = 0")

    @staticmethod
    def tanhResiduals(p, y, x):
        x0, d = p
        return y - 0.5 * (1 - numerix.tanh((x - x0) / (2*d)))

    def _plot(self):
        (x0_fit, d_fit), msg = leastsq(self.tanhResiduals,
                                       [self.L/2., self.deltaA],
                                       args=(self.phase.globalValue,
                                             self.x.globalValue))
        self.positions.append(x0_fit)
        self.times.append(float(self.elapsed))
        self.line.set_data(self.positions, self.times)
        self.bug.set_data([self.positions[-1]], [self.times[-1]])

        self.timer.set_text(r"$t = {:f}\,\mathrm{{s}}$".format(float(self.elapsed)))

        super(PhaseViewer, self)._plot()

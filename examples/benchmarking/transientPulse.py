#!/usr/bin/env python

if __name__ == "__main__":
    
    from benchmarker import Benchmarker
    bench = Benchmarker()

    bench.start()

    from fipy import *

    N = 100000
    L = 10.
    dx = L / N
    mesh = Grid1D(nx = N, dx = dx)

    bench.stop('mesh')

    bench.start()

    C = CellVariable(mesh = mesh)
    C.setValue(1, where=abs(mesh.getCellCenters()[0] - L/2.) < L / 10.)

    bench.stop('variables')

    bench.start()

    D = 1.

    eq = TransientTerm() == ImplicitDiffusionTerm(coeff = D)

    bench.stop('terms')

    ## from fipy import viewers
    ## viewer = viewers.make(vars = C, limits = {'datamin': 0, 'datamax': 1})
    ## viewer.plot()
    ## raw_input("initial")

    bench.start()

    dt = 1e0
    steps = 1
    for step in range(steps):
        eq.solve(var = C, dt = dt)
    ##     viewer.plot()

    bench.stop('solve')


    print bench.report(numberOfElements=N, steps=steps)

    ## raw_input("finished")

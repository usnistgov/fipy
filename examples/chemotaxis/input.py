"""

Input file for chemotaxis modeling.

Here are some test cases for the model.

>>> from __future__ import division

>>> from builtins import input
>>> from builtins import range
>>> from examples.chemotaxis.parameters import parameters
>>> from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm, ImplicitSourceTerm, Viewer, numerix

>>> params = parameters['case 2']

>>> nx = 50
>>> dx = 1.
>>> L = nx * dx

>>> mesh = Grid1D(nx=nx, dx=dx)

>>> shift = 1.

>>> KMVar = CellVariable(mesh=mesh, value=params['KM'] * shift, hasOld=1)
>>> KCVar = CellVariable(mesh=mesh, value=params['KC'] * shift, hasOld=1)
>>> TMVar = CellVariable(mesh=mesh, value=params['TM'] * shift, hasOld=1)
>>> TCVar = CellVariable(mesh=mesh, value=params['TC'] * shift, hasOld=1)
>>> P3Var = CellVariable(mesh=mesh, value=params['P3'] * shift, hasOld=1)
>>> P2Var = CellVariable(mesh=mesh, value=params['P2'] * shift, hasOld=1)
>>> RVar = CellVariable(mesh=mesh, value=params['R'], hasOld=1)

>>> PN = P3Var + P2Var

>>> KMscCoeff = params['chiK'] * (RVar + 1) * (1 - KCVar - KMVar.cellVolumeAverage)
>>> KMspCoeff = params['lambdaK'] / (1 + PN / params['kappaK'])
>>> KMEq = TransientTerm() - KMscCoeff + ImplicitSourceTerm(KMspCoeff)

>>> TMscCoeff = params['chiT'] * (1 - TCVar - TMVar.cellVolumeAverage)
>>> TMspCoeff = params['lambdaT'] * (KMVar + params['zetaT'])
>>> TMEq = TransientTerm() - TMscCoeff + ImplicitSourceTerm(TMspCoeff)

>>> TCscCoeff = params['lambdaT'] * (TMVar * KMVar).cellVolumeAverage
>>> TCspCoeff = params['lambdaTstar']
>>> TCEq = TransientTerm() - TCscCoeff + ImplicitSourceTerm(TCspCoeff)

>>> PIP2PITP = PN / (PN / params['kappam'] + PN.cellVolumeAverage / params['kappac'] + 1) + params['zetaPITP']

>>> P3spCoeff = params['lambda3'] * (TMVar + params['zeta3T'])
>>> P3scCoeff = params['chi3'] * KMVar * (PIP2PITP / (1 + KMVar / params['kappa3']) + params['zeta3PITP']) + params['zeta3']
>>> P3Eq = TransientTerm() - DiffusionTerm(params['diffusionCoeff']) - P3scCoeff + ImplicitSourceTerm(P3spCoeff)

>>> P2scCoeff = scCoeff = params['chi2'] + params['lambda3'] * params['zeta3T'] * P3Var
>>> P2spCoeff = params['lambda2'] * (TMVar + params['zeta2T'])
>>> P2Eq = TransientTerm() - DiffusionTerm(params['diffusionCoeff']) - P2scCoeff + ImplicitSourceTerm(P2spCoeff)

>>> KCscCoeff = params['alphaKstar'] * params['lambdaK'] * (KMVar / (1 + PN / params['kappaK'])).cellVolumeAverage
>>> KCspCoeff = params['lambdaKstar'] / (params['kappaKstar'] + KCVar)
>>> KCEq = TransientTerm() - KCscCoeff + ImplicitSourceTerm(KCspCoeff)

>>> eqs = ((KMVar, KMEq), (TMVar, TMEq), (TCVar, TCEq), (P3Var, P3Eq), (P2Var, P2Eq), (KCVar, KCEq))

>>> if __name__ == '__main__':
...     steps = 100
... else:
...     steps = 28

>>> for i in range(steps):
...     for var, eqn in eqs:
...         var.updateOld()
...     for var, eqn in eqs:
...         eqn.solve(var, dt=1.)

>>> accuracy = 1e-2
>>> print(KMVar.allclose(params['KM'], atol=accuracy))
1
>>> print(TMVar.allclose(params['TM'], atol=accuracy))
1
>>> print(TCVar.allclose(params['TC'], atol=accuracy))
1
>>> print(P2Var.allclose(params['P2'], atol=accuracy))
1
>>> print(P3Var.allclose(params['P3'], atol=accuracy))
1
>>> print(KCVar.allclose(params['KC'], atol=accuracy))
1

>>> PNView = PN / PN.cellVolumeAverage
>>> PNView.name = 'PN'

>>> KMView = KMVar / KMVar.cellVolumeAverage
>>> KMView.name = 'KM'

>>> TMView = TMVar / TMVar.cellVolumeAverage
>>> TMView.naem = 'TM'

>>> RVar[:] = params['S'] + (1 + params['S']) * params['G'] * numerix.cos((2 * numerix.pi * mesh.cellCenters[0]) / L)

>>> if __name__ == '__main__':
...     KMViewer = Viewer((PNView, KMView, TMView), title = 'Gradient Stimulus: Profile')
...
...     for i in range(100):
...         for var, eqn in eqs:
...             var.updateOld()
...         for var, eqn in eqs:
...             eqn.solve(var, dt=0.1)
...
...         KMViewer.plot()
...
...     input("finished")
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

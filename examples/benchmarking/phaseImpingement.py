#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "phaseImpingement.py"
 #                                     created: 1/18/06 {2:35:59 PM}
 #                                 last update: 7/5/07 {8:06:32 PM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 # Author: Daniel Wheeler
 # E-mail: daniel.wheeler@nist.gov
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards and
 # Technology by employees of the Federal Government in the course of their
 # official duties.  Pursuant to title 17 Section 105 of the United States
 # Code this document is not subject to copyright protection and is in the
 # public domain.  phaseImpingement.py is an experimental work.  NIST assumes
 # no responsibility whatsoever for its use by other parties, and makes no
 # guarantees, expressed or implied, about its quality, reliability, or any
 # other characteristic.  We would appreciate acknowledgement if the 
 # document is used.
 # 
 # This document can be redistributed and/or modified freely provided that 
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and 
 #  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 # 
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 1940-01-18 JEG 1.0 original
 # 
 # ########################################################################
 ##

r"""

This example benchmarks the speed and memory usage of solving the
2D phase-field problem with grain impingement. Run:
    
    $ python setup.py efficiency_test
"""
__docformat__ = 'restructuredtext'

if __name__ == "__main__":
    
    from fipy import *
    from fipy.tools.parser import parse

    from benchmarker import Benchmarker
    bench = Benchmarker()

    numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 400)

    bench.start()

    steps = 10
    nx = int(sqrt(numberOfElements))
    ny = nx
    Lx = 2.5 * nx / 100.
    dx = Lx / nx

    mesh = Grid2D(dx,dx,nx,nx)

    bench.stop('mesh')

    timeStepDuration = 0.02
    phaseTransientCoeff = 0.1
    thetaSmallValue = 1e-6
    beta = 1e5
    mu = 1e3
    thetaTransientCoeff = 0.01
    gamma= 1e3
    epsilon = 0.008
    s = 0.01
    alpha = 0.015
    temperature = 10.

    bench.start()

    phase = CellVariable(
        name = 'PhaseField',
        mesh = mesh,
        value = 0.
        )

    theta = ModularVariable(
        name = 'Theta',
        mesh = mesh,
        value = -pi + 0.0001,
        hasOld = 1
        )

    x, y = mesh.getCellCenters()
    for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
                             (Lx, 0., -2. * pi / 3.), 
                             (0., Lx, -2. * pi / 3. + 0.3), 
                             (Lx, Lx,  2. * pi / 3.)):
        segment = (x - a)**2 + (y - b)**2 < (Lx / 2.)**2
        phase.setValue(1., where=segment)
        theta.setValue(thetaValue, where=segment)

    bench.stop('variables')

    bench.start()



    def buildPhaseEquation(phase, theta):

        mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)
        thetaMag = theta.getOld().getGrad().getMag()
        implicitSource = mPhiVar * (phase - (mPhiVar < 0))
        implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag

        return TransientTerm(phaseTransientCoeff) == \
                  ExplicitDiffusionTerm(alpha**2) \
                  - ImplicitSourceTerm(implicitSource) \
                  + (mPhiVar > 0) * mPhiVar * phase

    phaseEq = buildPhaseEquation(phase, theta)
    def buildThetaEquation(phase, theta):

        phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
        phaseModSq = phaseMod * phaseMod
        expo = epsilon * beta * theta.getGrad().getMag()
        expo = (expo < 100.) * (expo - 100.) + 100.
        pFunc = 1. + exp(-expo) * (mu / epsilon - 1.)

        phaseFace = phase.getArithmeticFaceValue()
        phaseSq = phaseFace * phaseFace
        gradMag = theta.getFaceGrad().getMag()
        eps = 1. / gamma / 10.
        gradMag += (gradMag < eps) * eps
        IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
        diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)

        thetaGradDiff = theta.getFaceGrad() - theta.getFaceGradNoMod()
        sourceCoeff = (diffusionCoeff * thetaGradDiff).getDivergence()

        return TransientTerm(thetaTransientCoeff * phaseModSq * pFunc) == \
                   ImplicitDiffusionTerm(diffusionCoeff) \
                   + sourceCoeff

    thetaEq = buildThetaEquation(phase, theta)

    bench.stop('terms')

    theta.updateOld()
    phase.updateOld()
    thetaEq.solve(theta, dt = timeStepDuration)
    phaseEq.solve(phase, dt = timeStepDuration)

    bench.start()

    ##from profiler import Profiler
    ##from profiler import calibrate_profiler
    ##fudge = calibrate_profiler(10000)
    ##profile = Profiler('profile-HEAD-i686', fudge=fudge)

    for i in range(steps):
        theta.updateOld()
        phase.updateOld()
        thetaEq.solve(theta, dt = timeStepDuration)
        phaseEq.solve(phase, dt = timeStepDuration)

    ##profile.stop()

    bench.stop('solve')

    print bench.report(numberOfElements=numberOfElements, steps=steps)

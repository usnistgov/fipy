#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/6/04 {4:36:39 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""
Electrochemical Phase Field input file

Build a mesh, variable, and diffusion equation with fixed (zero) flux
boundary conditions at the top and bottom and fixed value boundary
conditions at the left and right.

Iterates a solution and plots the result with gist.

Iteration is profiled for performance.

The dimensionless governing equations for the electrochemical phase field 
problem consist of the phase field governing equation

.. raw:: latex

   \begin{align}
       \frac{\partial\Phase}{\partial\Time} 
       &= -\Mobility{\Phase}\left[
	   \Interpolate'\left(\Phase\right) 
	   \sum_{j=1}^{\Components} \Concentration{j} \Delta\Chemical{j}{\Standard}
	   + \DoubleWell'\left(\Phase\right)
	   \sum_{j=1}^{\Components} \Concentration{j} \Barrier{j}
       \right. \nonumber
       \\
       &\qquad\left.
	   \vphantom{\sum_{j=1}^{\Components} \Concentration{j} \Delta\Chemical{j}{\Standard}}
	   - \Gradient{\Phase}\nabla^2\Phase
	   - \frac{\Dielectric'(\Phase)}{2} 
	       \left(\nabla\Potential\right)^{2}
       \right]
       \label{eq:Evolution:Ideal:Phase} 
   \end{align}
   
Poisson's equation

.. raw:: latex

   \begin{equation}
       \nabla\cdot\left[\Dielectric(\Phase)\nabla\Potential\right] 
       + \sum_{j=1}^{\Components} \Valence{j} \Concentration{j} = 0
       \label{eq:Governing:Poisson}
   \end{equation}

and a diffusion equation for each species

.. raw:: latex

   \begin{equation}
       \frac{\partial\Concentration{j}}{\partial\Time}
       = - \nabla\cdot\Flux{j},
       \qquad j = 1 \ldots \Components-1
       \label{eq:Conservation:Concentration}
   \end{equation}

where the flux of substitional species is given by

.. raw:: latex

   \begin{equation}
       \Flux{j} 
       = -\Mobility{j}\nabla\left[	
	   \Delta\Chemical{j\Solvent{}}{\Standard}
		   \Interpolate\left(\Phase\right)
	   + \ln\frac{\Concentration{j}}{\Concentration{\Solvent{}}}
	   + \Valence{j\Solvent{}} \Potential
	   + \Barrier{j\Solvent{}} \DoubleWell\left(\Phase\right)
       \right]
       \qquad j = 2 \ldots \Components-1
       \label{eq:Flux:Ideal:Substitutional}
   \end{equation}

and the flux of electrons by

.. raw:: latex

   \begin{equation}
       \Flux{\Electron} 
	   = -\Mobility{\Electron}\nabla\left[
	       \Delta\Chemical{\Electron}{\Standard}\Interpolate\left(\Phase\right)
	       + \ln\frac{\Concentration{\Electron}}{1+\Concentration{\Electron}}
	       + \Valence{\Electron}\Potential
	       + \Barrier{\Electron} \DoubleWell\left(\Phase\right)
	   \right]
       \label{eq:Flux:Ideal:Interstitial}    
   \end{equation}
   
The mobilities are given by

.. raw:: latex

   \begin{equation}
       \Mobility{j} 
       = \frac{\Diffusivity{jj}\Concentration{\Solvent{}}\Concentration{j}}
	       {\Concentration{\Solvent{}}+\Concentration{j}}
	   \qquad j = 2 \ldots \Components-1
	   \label{eq:Mobility:Subsitutional}
   \end{equation}

and

.. raw:: latex

   \begin{equation}
       \Mobility{\Electron} 
       = \Diffusivity{\Electron}\left(
		   1+\Concentration{\Electron}
	       \right)\Concentration{\Electron}
	   \label{eq:Mobility:Electron}
   \end{equation}

We combine Eq.~\eqref{eq:Conservation:Concentration} with 
Eq.~\eqref{eq:Flux:Ideal:Substitutional} and rearranged to see the 
terms

.. raw:: latex

   \begin{equation}
       \underbrace{
	   \frac{\partial\Concentration{j}}{\partial\Time}
       }_\text{transient}
       = 
       \underbrace{
	   \nabla\cdot\Concentration{j}
	   \frac{\Diffusivity{jj}}{1-\sum_{\substack{i=2\\ i \neq j}}^{\Components-1} \Component{i}}
	   \left\{
	       \left(
		   1-\sum_{i=2}^{\Components-1} \Component{i}
	       \right)
	       \left(
		   \left[
		       \Delta\Chemical{j\Solvent{}}{\Standard}
			       \Interpolate'\left(\Phase\right)
		       + \Barrier{j\Solvent{}} \DoubleWell'\left(\Phase\right)
		   \right] \nabla\Phase
		   + \Valence{j\Solvent{}} \nabla\Potential
	       \right)
	       + \sum_{\substack{i=2\\ i \neq j}}^{\Components-1} \nabla\Component{i}
	   \right\}
       }_\text{convection}       
       + \underbrace{
	   \Diffusivity{jj}\nabla^2\Concentration{j}
       }_\text{diffusion}
   \end{equation}

"""

__docformat__ = 'restructuredtext'


import Numeric

## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D

from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.viewers.gist1DViewer import Gist1DViewer
from fipy.viewers.gistVectorViewer import GistVectorViewer
## from fipy.iterators.iterator import Iterator
## from fipy.iterators.adaptiveIterator import AdaptiveIterator
from fipy.models.elphf.myAdaptiveIterator import MyAdaptiveIterator

from fipy.variables.variable import Variable

import fipy.tools.dimensions.physicalField

import fipy.models.elphf.elphf as elphf

nx = 1000
dx = "0.00005 nm"
# L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = "1. nm",
    nx = nx,
    ny = 1)
    
parameters = {
    'time step duration': "1e-9 s",
    'substitutional molar volume': "1.80000006366754e-05 m**3/mol",
    'phase': {
	'name': "xi",
	'mobility': "1 m**3/J/s",
	'gradient energy': "3.6e-11 J/m",
	'value': 1.
    },
    'potential': {
	'name': "psi",
	'dielectric': 78.49
    },
    'solvent': {
	'standard potential': "34139.7265625 J/mol",
	'barrier height': "3.6e5 J/mol",
	'valence': 0
    }
}

parameters['interstitials'] = (
    {
	'name': "e-",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "-33225.9453125 J/mol",
	'barrier height': "0. J/mol",
	'valence': -1,
	'value': "111.110723815414 mol/l",
    },
)

parameters['substitutionals'] = (
    {
	'name': "SO4",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "24276.6640625 J/mol",
	'barrier height': parameters['solvent']['barrier height'],
	'valence': -2,
	'value': "0.000010414586295976 mol/l",
    },
    {
	'name': "Cu",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "-7231.81396484375 J/mol",
	'barrier height': parameters['solvent']['barrier height'],
	'valence': +2,
	'value': "55.5553718417909 mol/l",
    }
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(lambda cell: cell.getCenter()[0] > mesh.getPhysicalShape()[0]/2.)
fields['phase'].setValue(0.,setCells)
fields['interstitials'][0].setValue("0.000111111503177394 mol/l", setCells)
fields['substitutionals'][0].setValue("0.249944439430068 mol/l", setCells)
fields['substitutionals'][1].setValue("0.249999982581341 mol/l", setCells)

## fields['substitutionals'][0][nx/2-5] = 0.6
## fields['phase'][nx/2-9] = 0.9
## fields['phase'][nx/2-8] = 0.8
## fields['phase'][nx/2-7] = 0.7
## fields['phase'][nx/2-6] = 0.6
## fields['phase'][nx/2-5] = 0.5
## fields['phase'][nx/2-4] = 0.4
## fields['phase'][nx/2-3] = 0.3
## fields['phase'][nx/2-2] = 0.2
## fields['phase'][nx/2-1] = 0.1
## fields['substitutionals'][1][nx/2-5] = 1.5


phaseViewer = Gist1DViewer(vars = (fields['phase'],))
## phaseViewer = Gist1DViewer(vars = (fields['phase'].get_gPrime(),))
## phaseViewer1 = GistVectorViewer(var = fields['phase'].get_p().getFaceGrad() )
## phaseViewer2 = GistVectorViewer(var = fields['phase'].get_pPrime().getArithmeticFaceValue().transpose()*fields['phase'].getFaceGrad())
## phaseViewer3 = Gist1DViewer(vars = (fields['phase'].get_pPrime().getArithmeticFaceValue(),))
potentialViewer = Gist1DViewer(vars = (fields['potential'],))
## concViewer = Gist1DViewer(vars = list(fields['substitutionals']) + list(fields['interstitials']) + [fields['solvent']], ylog = 1)
concViewer = Gist1DViewer(vars = list(fields['substitutionals']) + list(fields['interstitials']) + [fields['solvent']])
## concViewer = Gist1DViewer(vars = [field.getGrad() for field in list(fields['substitutionals'])] + [field.getGrad for field in list(fields['interstitials'])] + [fields['solvent'].getGrad()])

## phaseRelaxation = Variable(value = 1.)
phaseRelaxation = 1.

equations = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters,
    phaseRelaxation = phaseRelaxation
)

## timeStepDuration = fipy.tools.dimensions.physicalField.Scale(parameters['time step duration'], "TIME")

chargeViewer = Gist1DViewer(vars = (fields['charge'],))

viewers = (phaseViewer, potentialViewer, concViewer,chargeViewer)
	
## it = Iterator(equations = equations)
it = MyAdaptiveIterator(equations = equations, viewers = viewers, phaseRelaxation = phaseRelaxation)

if __name__ == '__main__':
    for viewer in viewers:
	viewer.plot()

    print fields['substitutionals'][1].name, Numeric.sum(fields['substitutionals'][1].getNumericValue()), fields['substitutionals'][1][540:560]
    
    raw_input()

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    # it.timestep(50)
    # 
    # for timestep in range(5):
    #     it.sweep(50)
    #     it.advanceTimeStep


    for i in range(50):
	try:
	    it.timestep(maxSweeps = 20, dt = 1.)
## 	    it.timestep(dt = 1.)

            for viewer in viewers:
                viewer.plot()
		
	    print fields['substitutionals'][1].name, Numeric.sum(fields['substitutionals'][1].getNumericValue()), fields['substitutionals'][1][540:560]
		
    ## 	it.timestep(steps = 1, maxSweeps = 5)
	except KeyboardInterrupt:
	    break
	except Exception, e:
            raise
	    print "Error:", e
	except:
	    print "Not converged"
	    
	print "***** timestep", i, "******"

        raw_input()

	
    raw_input()
	
    # profile.stop()
	    
    ## raw_input()


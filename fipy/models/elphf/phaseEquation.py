#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "phaseEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 2/18/05 {10:41:27 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

r"""
Represents the phase field equation

.. raw:: latex

   \[
       \underbrace{
	   \frac{\partial \xi}{\partial t}
	   \vphantom{M_{\xi}\sum_{j=1}^{n} C_j}
       }_{\text{transient}}
       = 
       \underbrace{
	   M_{\xi}\kappa_{\xi}\nabla^2 \xi
	   \vphantom{M_{\xi}\sum_{j=1}^{n} C_j}
       }_{\text{diffusion}}
       - 
       \underbrace{
	   \overbrace{
	       M_{\xi}\sum_{j=1}^{n} C_j \left[
		   p'(\xi) \Delta\mu_j^\circ
		   + g'(\xi) W_j
	       \right]
	   }^{\text{phase transformation}}
	   % +
	   % \overbrace{
	   %   M_{\xi}\frac{\epsilon'(\xi)}{2}\left(\nabla\phi\right)^2
	   % }^{\text{dielectric}}
       }_{\text{source}}
   \]

   where \( \xi \) is the phase field variable, \( t \) is time, \(
   M_\xi \) is the phase field mobility, \( \kappa_\xi \) is the phase
   field gradient energy coefficient.  \( p(\xi) \) describes the
   interpolation of the free energy, \( g(\xi) \) describes the
   shape of the energy barrier between the electrode and electrolyte
   phases, \( p'(\xi) = 30\xi^2\left(1-\xi\right)^2 \), and \( g'(\xi) =
   2\xi\left(1-\xi\right)\left(1-2\xi\right) \).  For a given species
   \( j \), \( C_j \) is the concentration, \( \Delta\mu_j^{\circ} \)
   is the standard chemical potential difference between the electrode
   and electrolyte for a pure material, and \( W_j \) is the magnitude
   of the energy barrier in the double-well free energy function.
"""
__docformat__ = 'restructuredtext'

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.tools.dimensions import physicalField

from equationFactory import EquationFactory

class PhaseEquationFactory(EquationFactory):
    def make(self, fields, parameters):
	from elphf import constant

	phaseMobility = physicalField.Scale(parameters['phase']['mobility'], 
					    constant['MOLARVOLUME'] 
					    / (constant['ENERGY'] * constant['TIME']))
	phaseGradientEnergy = physicalField.Scale(parameters['phase']['gradient energy'], 
						  constant['LENGTH']**2 * constant['ENERGY'] 
						  / constant['MOLARVOLUME'])
						  
	enthalpy = fields['solvent'].getStandardPotential()
	barrier = fields['solvent'].getBarrierHeight()

	for component in list(fields['substitutionals']) + list(fields['interstitials']):
	    enthalpy = enthalpy + component * component.getStandardPotential() #.getOld()
	    barrier = barrier + component * component.getBarrierHeight() #.getOld()

	mPhi = -(30. * fields['phase'] * (1. - fields['phase']) * enthalpy \
	    + 2 * (1. - 2 * fields['phase']) * barrier)

	return TransientTerm(tranCoeff = 1. / phaseMobility) \
	    - ImplicitDiffusionTerm(diffCoeff = phaseGradientEnergy) \
	    + ImplicitSourceTerm(sourceCoeff = mPhi * (fields['phase'] - (mPhi < 0.))) \
            - (mPhi > 0.) * mPhi * fields['phase']

factory = PhaseEquationFactory()
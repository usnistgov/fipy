#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "interstitialEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 4/1/05 {9:28:09 PM} 
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
Represents the diffusion equation for interstitial species, such as 
electrons,

.. raw:: latex

   \[
       \underbrace{
	   \frac{\partial C_j}{\partial t}
	   \vphantom{\left\{
	       \overbrace{
		   \left[p'(\xi)\right]
	       }^{\text{phase transformation}}
	   \right\}}
       }_{\text{transient}}
       = \underbrace{
	   D_j\nabla^2 C_j
	   \vphantom{\left\{
	       \overbrace{
		   \left[p'(\xi)\right]
	       }^{\text{phase transformation}}
	   \right\}}
       }_{\text{diffusion}} \\
       + \underbrace{
	   D_j\nabla\cdot 
	   C_j
	   \left\{
	       \overbrace{
		   \left[
		       p'(\xi) \Delta\mu_j^{\circ}
		       + g'(\xi) W_j
		   \right] \nabla\xi
	       }^{\text{phase transformation}}
	       +
	       \overbrace{
		   z_j \nabla \phi
	       }^{\text{electromigration}}
	   \right\}
       }_{\text{convection}}
   \]

   where, for a given species \( j \), \( C_j \) is the concentration, 
   \( D_j \) is the self diffusivity, \( \Delta\mu_j^{\circ} \) is the 
   standard chemical potential difference between the electrode and 
   electrolyte for a pure material, \( W_j \) is the magnitude of 
   the energy barrier in the double-well free energy function, and \( z_j \)
   is the valence.
   
   In addition, \( t \) is time, \( \xi \) is the phase field variable, 
   \( \phi \) is the electrostatic potential,
   \( p(\xi) \) describes the
   interpolation of the free energy, \( g(\xi) \) describes the
   shape of the energy barrier between the electrode and electrolyte
   phases, \( p'(\xi) = 30\xi^2\left(1-\xi\right)^2 \), and \( g'(\xi) =
   2\xi\left(1-\xi\right)\left(1-2\xi\right) \).  

"""
__docformat__ = 'restructuredtext'

from concentrationEquation import ConcentrationEquationFactory

class InterstitialEquationFactory(ConcentrationEquationFactory):
    def getConvectionCoeff(self, Cj, fields, diffusivity = None):
	if diffusivity is None:
	    diffusivity = Cj.getDiffusivity()
	Cj.weightedDiffusivity = (diffusivity * (1. + Cj.getHarmonicFaceValue())).transpose()
	
	return ConcentrationEquationFactory.getConvectionCoeff(self, Cj = Cj, fields = fields, 
							       diffusivity = Cj.weightedDiffusivity)
							       
factory = InterstitialEquationFactory()
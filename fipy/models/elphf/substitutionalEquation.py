#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "substitutionalEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 12/8/04 {5:10:20 PM} 
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
Represents the diffusion equation for substitutional species

.. raw:: latex

   \begin{align*}
       \underbrace{
	   \frac{\partial C_j}{\partial t}
       }_{\text{transient}}
       &= \underbrace{
	   D_{j}\nabla^2 C_j
	   \vphantom{\frac{\partial C_j}{\partial t}}
       }_{\text{diffusion}} \\
       & \qquad + \underbrace{
	   D_{j}\nabla\cdot 
	   \frac{C_j}{1 - \sum_{\substack{k=2\\ k \neq j}}^{n-1} C_k}
	   \left\{
	       \overbrace{
		   \sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i
	       }^{\text{counter diffusion}}
	       + 
	       \overbrace{
		   C_n \left[
		       p'(\xi) \Delta\mu_{jn}^{\circ}
		       + g'(\xi) W_{jn}
		   \right] \nabla\xi
		   \vphantom{\sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i}
	       }^{\text{phase transformation}}
	       +
	       \overbrace{
		   C_n z_{jn} \nabla \phi
		   \vphantom{\sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i}
	       }^{\text{electromigration}}
	   \right\}
       }_{\text{convection}}
   \end{align*}
   
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
   
   The summation \( \sum_{\substack{i=2\\ i \neq j}}^{n-1} \) is over 
   all substitutional species, excluding the species of interest and the 
   designated solvent, and \( \Delta\mu_{jn}^{\circ} \), \( W_{jn} \), and 
   \( z_{jn} \) are the differences of the respective quantities 
   \( \Delta\mu_{j}^{\circ} \), \( W_{j} \), and \( z_{j} \) between
   substitutional species \( j \) and the solvent species \( n \).
"""
__docformat__ = 'restructuredtext'

from concentrationEquation import ConcentrationEquationFactory

class SubstitutionalEquationFactory(ConcentrationEquationFactory):
    def getConvectionCoeff(self, Cj, fields, diffusivity = None):
	Cj.substitutionalSum = Cj.copy()
        Cj.substitutionalSum.setValue(0)

	for component in [component for component in fields['substitutionals'] if component is not Cj]:
	    Cj.substitutionalSum = Cj.substitutionalSum + component

	denom = 1. - Cj.substitutionalSum.getHarmonicFaceValue()
	if diffusivity is None:
	    diffusivity = Cj.getDiffusivity()
	    
	Cj.subsConvCoeff = diffusivity * Cj.substitutionalSum.getFaceGrad() / denom.transpose()
	Cj.weightedDiffusivity = (diffusivity * fields['solvent'].getHarmonicFaceValue() / denom).transpose()

	return Cj.subsConvCoeff \
	    + ConcentrationEquationFactory.getConvectionCoeff(self, Cj = Cj, fields = fields, 
							      diffusivity = Cj.weightedDiffusivity)
	
factory = SubstitutionalEquationFactory()
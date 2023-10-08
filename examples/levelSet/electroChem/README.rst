=========================================
Superconformal Electrodeposition Examples
=========================================

The Damascene Process
---------------------

State of the art manufacturing of semiconductor devices involves the
electrodeposition of copper for on-chip wiring of integrated circuits.
In the Damascene process interconnects are fabricated by first
patterning trenches in a dielectric medium and then filling by metal
electrodeposition over the entire wafer surface. This metallization
process, pioneered by IBM, depends on the use of electrolyte additives
that effect the local metal deposition rate.

Superfill
---------

The additives in the electrolyte affect the local deposition rate in
such a way that bottom-up filling occurs in trenches or vias. This
process, known as superconformal electrodeposition or superfill, is
demonstrated in the following figure. The figure shows sequential
images of bottom-up superfilling of submicrometer trenches by copper
deposition from an electrolyte containing PEG-SPS-Cl. Preferential
metal deposition at the bottom of the trenches followed by bump
formation above the filled trenches is evident.

.. image:: /figures/examples/levelSet/electroChem/experimentalSuperfill.*
   :scale: 50
   :align: center
   :alt: cross-sectional scanning electron micrographs of experimental superfill

The CEAC Mechanism
------------------

This process has been demonstrated to depend critically on the
inclusion of additives in the electrolyte.  Recent publications
propose Curvature Enhanced Accelerator Coverage (CEAC) as the
mechanism behind the superfilling process :cite:`NIST:damascene:2001`.  In this
mechanism, molecules that accelerate local metal deposition displace
molecules that inhibit local metal deposition on the metal/electrolyte
interface. For electrolytes that yield superconformal filling of fine
features, this buildup happens relatively slowly because the
concentration of accelerator species is much more dilute compared to
the inhibitor species in the electrolyte.  The mechanism that leads to
the increased rate of metal deposition along the bottom of the filling
trench is the concurrent local increase of the accelerator coverage
due to decreasing local surface area, which scales with the local
curvature (hence the name of the mechanism). A good overview of this
mechanism can be found in :cite:`moffatInterface:2004`.

Using FiPy to model Superfill
-----------------------------

Example :mod:`~examples.levelSet.electroChem.simpleTrenchSystem` provides a simple way to use :term:`FiPy` to model the
superfill process. The example includes a detailed description of the governing
equations and feature geometry. It requires the user to import and execute a
function at the python prompt.  The model parameters can be passed as arguments
to this function. In future all superfill examples will be provided with this
type of interface. Example :mod:`~examples.levelSet.electroChem.howToWriteAScript` has the same functionality as
:mod:`~examples.levelSet.electroChem.simpleTrenchSystem` but demonstrates how to write a new script in the case where
the existing suite of scripts do not meet the required needs.

In general it is a good idea to obtain the Mayavi_ plotting package
for which a specialized superfill viewer class has been created, see
:ref:`INSTALLATION`. The other standard viewers mentioned in
:ref:`INSTALLATION` are still adequate although they do not give such
clear images that are tailored for the superfill problem.  The images
below demonstrate the Mayavi_ viewing capability.  Each contour
represents sequential positions of the interface and the color
represents the concentration of accelerator as a surfactant. The areas
of high surfactant concentration have an increased deposition rate.

.. image:: /figures/examples/levelSet/electroChem/superfillImage.*
   :scale: 50
   :align: left
   :alt: FiPy simulation contours of superfill

.. _Mayavi:   http://mayavi.sourceforge.net

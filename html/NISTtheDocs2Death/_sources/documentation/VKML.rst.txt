========================================
Virtual Kinetics of Materials Laboratory
========================================

The VKML is a set of simple :term:`FiPy` examples that
simulate basic aspects of kinetics of materials through an interactive
Graphical User Interface. The seminal development by Michael Waters and
`Prof. R. Edwin Garcia`_ of `Purdue University`_ includes four examples:

 `Polycrystalline Growth and Coarsening`_

   simulates the growth, impingement, and coarsening of a random
   distribution of crystallographically oriented nuclei. The user can
   control every aspect of the model such as the nuclei radius, the size of
   the simulation cell, and whether the grains are homogeneously dispersed
   or only on one wall of the simulation.

 `Dendritic Growth`_

   simulates the anisotropic solidification of a single solid seed with an
   N-fold axis of crystallographic symmetry embedded in an undercooled
   liquid. The user can specify many material aspects of the solidification
   process, such as the thermal diffusivity and the strength of the surface
   tension anisotropy. Default values are physical but arbitrary. This
   model is based on the phase field method and an example shown in the
   :term:`FiPy` manual.

 `Two-Dimensional Spinodal Decomposition`_

   simulates the time-dependent segregation of two chemical components and
   its subsequent coarsening, as presented by John Cahn. The default values
   are physical but arbitrary.

 `Three-Dimensional Spinodal Decomposition`_

   has the same functionality as the 2D version, but has an interactive
   Three-Dimensional viewer.

These modules provide a Graphical User Interface to :term:`FiPy`, and allow you to
perform simulations directly through your web browser. This approach to
computing removes the need to install the software on your local machine
(unless you really want to), allows you to assess current and potential
:term:`FiPy` applications and instead you only need a web browser to access it and
run it. In other words, you can run these simulations (and simulations like
this one) from a Windows machine, a Mac, or a Linux box, and you can also
run the modules from Michigan, Boston, Japan, or England: from wherever you
are. Moreover, if you close your web browser and leave your calculation
running, when you come back a few hours later, your calculation will
persist. Additionally, if there is something you want to share with a
coworker, wherever he or she might be (e.g., the other side of the planet),
you can grant him temporary access to your calculation so that the third
party can directly see the output (or specify inputs directly into it,
without having to travel to where you are). It is a great way to privately
(or publicly) collaborate with other people, especially if the users are in
different parts of the world.

The only requirement to run VKML is to register (registration is 100% free)
in the nanoHUB_.


.. _Prof. R. Edwin Garcia:                    http://bicephalous.ecn.purdue.edu/~edwin
.. _Purdue University:                        http://www.purdue.edu
.. _Polycrystalline Growth and Coarsening:    https://www.nanohub.org/tools/vkmlpsgg/
.. _Dendritic Growth:                         https://www.nanohub.org/tools/vkmlggs/
.. _Two-Dimensional Spinodal Decomposition:   https://www.nanohub.org/tools/vkmlsd/
.. _Three-Dimensional Spinodal Decomposition: https://www.nanohub.org/tools/vkmlsd3d/

.. _nanoHUB:                                  http://www.nanoHUB.org

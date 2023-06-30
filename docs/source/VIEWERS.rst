.. _VIEWERS:

=======
Viewers
=======

A viewer is required to see the results of :term:`FiPy`
calculations. :term:`Matplotlib` is by far the most widely used
:term:`Python` based viewer and the best choice to get :term:`FiPy` up
and running quickly. :term:`Matplotlib` is also capable of publication
quality plots. :term:`Matplotlib` has only rudimentary 3D capability,
which :term:`FiPy` does not attempt to use. :term:`Mayavi` is required for
3D viewing.

.. _MATPLOTLIB:

----------
Matplotlib
----------

http://matplotlib.sourceforge.net

:term:`Matplotlib` is a :term:`Python` package that displays
publication quality results. It displays both 1D X-Y type plots and 2D
contour plots for both structured and unstructured data, but does not
display 3D data. It works on all common platforms.

.. _MAYAVI:

------
Mayavi
------

http://code.enthought.com/projects/mayavi/

The :term:`Mayavi` Data Visualizer is a free, easy to use scientific data
visualizer.  It displays 1D, 2D and 3D data. It is the only
:term:`FiPy` viewer available for 3D data. :ref:`MATPLOTLIB` is
probably a better choice for 1D or 2D viewing.

:term:`Mayavi` requires VTK_, which can be difficult to build from source.

.. note::

   MayaVi 1 is no longer supported.

.. _VTK: http://www.vtk.org/
.. _Mac OS X: http://www.apple.com/macosx

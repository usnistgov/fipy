#!/usr/bin/env python
# coding: utf-8

# # Phase Field Benchmark 8c
# ## Explicit nucleation, multiple seeds at random times
# FiPy implementation of problem 3 in *Nucleation Benchmark Problem*, Wenkun Wu *et al.*, January 16, 2020
# 
# Based on problem 2.3 in *Benchmark problems for nucleation*, Tamás Pusztai, September 25, 2019

# **Do not edit `benchmark8c.py`**. Generate the batch-runnable file from the notebook with
# ```bash
# jupyter nbconvert benchmark8c.ipynb --to python
# ```

# ## Import Python modules

# In[ ]:


import os
import re
import sys
import yaml

import datreant as dtr

import fipy as fp
from fipy.tools import parallelComm
from fipy.meshes.factoryMeshes import _dnl


# Jupyter notebook handles some things differently than from the commandline

# In[ ]:


try:
    from IPython import get_ipython
    isnotebook = get_ipython() is not None
except:
    isnotebook = False


# ## Initialize
# ### Load parameters

# In[ ]:


if isnotebook:
    yamlfile = "params8c.yaml"
else:
    yamlfile = sys.argv[1]

with open(yamlfile, 'r') as f:
    if hasattr(yaml, "FullLoader"):
        # PyYAML 5.1 deprecated the plain yaml.load(input) function
        # https://github.com/yaml/pyyaml/wiki/PyYAML-yaml.load(input)-Deprecation
        params = yaml.load(f, Loader=yaml.FullLoader)
    else:
        params = yaml.load(f)


# ### Set any parameters for interactive notebook

# In[ ]:


if isnotebook:
    params['Lx'] = 100.
    params['Ly'] = 100.
    params['checkpoint_interval'] = 1.5 * params['dt']
    params['savetime'] = 4 * params['dt'] 
    params['totaltime'] = 100 * params['dt']
#    params['restart'] = "t=0.0.tar.gz"


# ### Initialize mesh and solution variables
# 
# Either restart from some `path/to/t={time}.tar.gz`, where the time is assigned to `elapsed`
# 
# or
# 
# Create a mesh based on parameters. Set
# >  the domain size to 1000 × 1000... the spatial and temporal resolution by setting $\Delta x = \Delta y = 0.8$ and $\Delta t = 0.04$

# In[ ]:


checkpoint_interval = params['checkpoint_interval']
savetime = params['savetime']
totaltime = params['totaltime']
dt = params['dt']

if params['restart']:
    phi, = fp.tools.dump.read(filename=params['restart'])
    mesh = phi.mesh

    X, Y = mesh.faceCenters
    
    Lx = mesh.communicator.MaxAll(max(X)) - mesh.communicator.MinAll(min(X))
    Ly = mesh.communicator.MaxAll(max(Y)) - mesh.communicator.MinAll(min(Y))

    # scanf("%g") simulator
    # https://docs.python.org/3/library/re.html#simulating-scanf
    scanf_g = "[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    pattern = ".*t=({g})\.tar\.gz".format(g=scanf_g)
    elapsed = re.match(pattern, params['restart']).group(1)
    elapsed = fp.Variable(name="$t$", value=float(elapsed))
else:
    Lx = params['Lx']
    Ly = params['Ly']

    dx, nx = _dnl(dx=params['dx'], nx=None, Lx=Lx)
    dy, ny = _dnl(dx=params['dx'], nx=None, Lx=Ly)

    mesh = fp.Grid2D(dx=dx, nx=nx, dy=dy, ny=ny)

    phi = fp.CellVariable(mesh=mesh, name="$\phi$", value=0., hasOld=True)

    elapsed = fp.Variable(name="$t$", value=0.)
    
x, y = mesh.cellCenters[0], mesh.cellCenters[1]
X, Y = mesh.faceCenters[0], mesh.faceCenters[1]


# In[ ]:


if isnotebook:
    viewer = fp.Viewer(vars=phi, datamin=0., datamax=1.)
    viewer.plot()


# ## Define governing equation

# > use only the nondimensional forms of the phase-field and nucleation equations, but without the tildes, for simplicity
# 
# > [Set] the driving force to $\Delta f = 1 / (6\sqrt{2})$

# In[ ]:


Delta_f = 1. / (6 * fp.numerix.sqrt(2.))


# > $$r_c = \frac{1}{3\sqrt{2}}\frac{1}{\Delta f} = 2.0$$

# In[ ]:


rc = 2.0


# and define the governing equation 
# > \begin{align}
# \frac{\partial\phi}{\partial t} &= \nabla^2\phi - g'(\phi) + \Delta f p'(\phi) \tag{7}
# \end{align}
# 
# > $$g(\phi) = \phi^2(1 - \phi)^2 \qquad p(\phi)=\phi^3(10 - 15\phi + 6\phi^2)$$
# 
# Following [`examples/phase/simple.py`](https://www.ctcms.nist.gov/fipy/examples/phase/generated/examples.phase.simple.html)
# 
# 
# \begin{align}
# \frac{\partial\phi}{\partial t}
# &= \nabla^2\phi + m_\phi \phi (1 - \phi) \notag
# \qquad\text{for $m_\phi \equiv -2(1 - 2\phi) + 30 \phi (1 - \phi) \Delta f$} \notag
# \\
# &= \nabla^2\phi + S \notag
# \\
# &\approx \nabla^2\phi + S|_\text{old}
# + \left.{\frac{\partial S}{\partial \phi}}\right|_\text{old} (\phi - \phi_\text{old}) 
# \notag \\
# &= \nabla^2\phi + \left(S - \frac{\partial S}{\partial \phi} \phi\right)_\text{old} 
# + \left.{\frac{\partial S}{\partial \phi}}\right|_\text{old} \phi \notag
# \\
# &= \nabla^2\phi + S_0 + S_1 \phi \notag
# \\
# S_0 &\equiv \left(S - \frac{\partial S}{\partial \phi} \phi\right)_\text{old}
# \notag \\
# &= m_\phi \phi_\text{old} (1 - \phi_\text{old}) - S_1 \phi_\text{old}
# \notag \\
# S_1 &\equiv \left.{\frac{\partial S}{\partial \phi}}\right|_\text{old}
# \notag \\
# &= \frac{\partial m_\phi}{\partial \phi} \phi (1 - \phi) + m_\phi (1 - 2\phi)
# \notag
# \end{align}

# In[ ]:


mPhi = -2 * (1 - 2 * phi) + 30 * phi * (1 - phi) * Delta_f
dmPhidPhi = 4 + 30 * (1 - 2 * phi) * Delta_f
S1 = dmPhidPhi * phi * (1 - phi) + mPhi * (1 - 2 * phi)
S0 = mPhi * phi * (1 - phi) - S1 * phi

eq = (fp.TransientTerm() == 
      fp.DiffusionTerm(coeff=1.) + S0 + fp.ImplicitSourceTerm(coeff=S1))


# ## Calculate total free energy
# 
# > \begin{align}
# F[\phi] = \int\left[\frac{1}{2}(\nabla\phi)^2 + g(\phi) - \Delta f p(\phi)\right]\,dV \tag{6}
# \end{align}

# In[ ]:


ftot = (0.5 * phi.grad.mag**2
        + phi**2 * (1 - phi)**2
        - Delta_f * phi**3 * (10 - 15 * phi + 6 * phi**2))
volumes = fp.CellVariable(mesh=mesh, value=mesh.cellVolumes)
F = ftot.cellVolumeAverage * volumes.sum()


# ## Define nucleation
# 
# > generate ... supercritical seeds with $r_0 = 1.1r^∗$. ... When adding a new seed, simply add the $\phi$ values given by the $\phi(r)$ profile
# \begin{align}
# \phi(x) &= \frac{1 - \tanh\left(\frac{x - x_0}{\sqrt{2}}\right)}{2}\tag{8}
# \end{align}
# to the $\phi$ values already in the domain, and handle the possible overlaps by setting $\phi = 1$ for all cells where $\phi > 1.$

# In[ ]:


def nucleus(x0, y0, r0):
    r = fp.numerix.sqrt((x - x0)**2 + (y - y0)**2)

    return (1 - fp.numerix.tanh((r - r0) / fp.numerix.sqrt(2.))) / 2


# ### Determine nucleation times
# Either load nucleation times from `path/to/restart/nucleii.txt`, based on directory of `params['restart']`
# 
# or
# 
# > generate 100 random nucleation times in the range $t=0\dots100$ for adding the 100 seeds to the simulation domain.

# In[ ]:


if parallelComm.procID == 0:
    if params['restart']:
        fname = os.path.join(os.path.dirname(params['restart']), "nucleii.txt")
        nucleii = fp.numerix.loadtxt(fnamem, skiprows=1)
    else:
        times = fp.numerix.random.random(params['numnuclei']) * totaltime
        times.sort()
        nucleii = fp.numerix.concatenate((times[..., fp.numerix.newaxis],
                                          fp.numerix.random.random((params['numnuclei'], 2))),
                                         axis=-1)
else:
    nucleii = None
nucleii = parallelComm.bcast(nucleii, root=0)


# ## Setup output

# ### Setup ouput storage

# In[ ]:


try:
    from sumatra.projects import load_project
    project = load_project(os.getcwd())
    record = project.get_record(params["sumatra_label"])
    output = record.datastore.root
except:
    # either there's no sumatra, no sumatra project, or no sumatra_label
    # this will be the case if this script is run directly
    output = os.getcwd()
    
if parallelComm.procID == 0:
    print("storing results in {0}".format(output))
    data = dtr.Treant(output)
else:
    class dummyTreant(object):
        categories = dict()

    data = dummyTreant()


# ### Create particle counter

# In[ ]:


from scipy import ndimage

class LabelVariable(fp.CellVariable):
    """Label features in `var` using scipy.ndimage.label
    
    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        Field to be labeled. Any values equal to or greater than `threshold`
        are counted as features and values below are considered the background.
        
        .. important:
           Only sensible if `var` is defined on a `...Grid...` Mesh.
    structure : array_like, optional
        A structuring element that defines feature connections.
        `structure` must be centrosymmetric
        (see ```scipy.ndimage.label`` Notes
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html#scipy.ndimage.label>`_).
        If no structuring element is provided,
        one is automatically generated with a squared connectivity equal to
        one.  That is, for a 2-D `input` array, the default structuring element
        is::
            [[0,1,0],
             [1,1,1],
             [0,1,0]]
    threshold : float, optional
        Boundary value between features (inclusive) and background.
    dtype : date-type, optional
        The desired data-type for the labels. Note that the type must be able
        to store the largest label, or this Variable will raise an Exception.
        Default: int.
    """
    def __init__(self, var, structure=None, threshold=0.5, dtype=int):
        # We want our value to hold dtype,
        # but if we pass an array, the CellVariable
        # will probably be wonky
        value = fp.numerix.array(0.).astype(dtype).item()
        fp.CellVariable.__init__(self,
                                 mesh=var.mesh,
                                 value=value,
                                 elementshape=var.shape[:-1])
        self.var = self._requires(var)
        self.structure = structure
        self.threshold = threshold
        self.dtype = dtype
        self._num_features = None
    
    def _calcValue(self):
        """Label features of `var`
        
        Side-effect: sets self._num_features
        """
        arr = (self.var.globalValue > self.threshold).astype(self.dtype)
        shape = (self.var.mesh.args['nx'], self.var.mesh.args['ny'])
        arr = arr.reshape(shape)
        self._num_features = ndimage.label(input=arr,
                                           structure=self.structure,
                                           output=arr)
        return arr.flat
        
    @property
    def num_features(self):
        """How many objects were found
        """
        if self.stale or not self._isCached() or self._num_features is None:
            self._getValue()

        return self._num_features


# In[ ]:


labels = LabelVariable(phi, threshold=0.5)


# In[ ]:


if isnotebook:
    labelViewer = fp.Viewer(vars=labels, datamin=0, datamax=25)
    labelViewer.plot()


# ### Define output routines

# In[ ]:


def saveStats(elapsed):
    if parallelComm.procID == 0:
        fname = data['stats.txt'].make().abspath
        if os.path.exists(fname):
            # backup before overwrite
            os.rename(fname, fname + ".save")
        try:
            fp.numerix.savetxt(fname,
                               stats,
                               delimiter="\t",
                               comments='',
                               header="\t".join(["time", "fraction", "particle_count", "energy"]))
        except:
            # restore from backup
            os.rename(fname + ".save", fname)
        if os.path.exists(fname + ".save"):
            os.remove(fname + ".save")

def current_stats(elapsed):
    return [float(x) for x in [elapsed, phi.cellVolumeAverage, labels.num_features, F]]

def savePhi(elapsed):
    if parallelComm.procID == 0:
        fname = data["t={}.tar.gz".format(elapsed)].make().abspath
    else:
        fname = None
    fname = parallelComm.bcast(fname)

    fp.tools.dump.write((phi,), filename=fname)

def checkpoint(elapsed):
    saveStats(elapsed)
    savePhi(elapsed)


# ### Figure out when to save

# In[ ]:


checkpoints = (fp.numerix.arange(int(elapsed / checkpoint_interval),
                                 int(totaltime / checkpoint_interval)) + 1) * checkpoint_interval
for sometime in [savetime, totaltime]:
    if sometime > elapsed and sometime not in checkpoints: 
        checkpoints = fp.tools.concatenate([checkpoints, [sometime]])
checkpoints.sort()


# ### Output initial condition

# In[ ]:


if params['restart']:
    fname = os.path.join(os.path.dirname(params['restart']), "stats.txt")
    stats = fp.numerix.loadtxt(fname, skiprows=1)
    stats = stats[stats[..., 0] <= elapsed].tolist()
else:
    stats = []
    stats.append(current_stats(elapsed))

checkpoint(elapsed)
    
if parallelComm.procID == 0:
    fp.numerix.savetxt(data['nucleii.txt'].make().abspath, nucleii, 
                       delimiter="\t", comments='',
                       header="\t".join(["time", "x", "y"]))


# ## Solve and output

# In[ ]:


times = fp.tools.concatenate([checkpoints, nucleii[..., 0]])
times.sort()
times = times[(times > elapsed) & (times <= totaltime)]


# In[ ]:


for until in times:
    while elapsed.value < until:
        phi.updateOld()
        dt_until = (until - elapsed).value
        dt_save = dt
        if dt_until < dt:
            dt = dt_until
        for sweep in range(5):
            eq.sweep(var=phi, dt=dt)
        elapsed.value = elapsed() + dt
        stats.append(current_stats(elapsed))
        dt = dt_save

    for tt, fx, fy in nucleii[nucleii[..., 0] == until]:
        phi.setValue(phi + nucleus(x0=fx * Lx, y0=fy * Ly, r0=params['factor'] * 2))
        phi.setValue(1., where=phi > 1.)
              
    if elapsed in checkpoints:
        checkpoint(elapsed)
              
    if isnotebook:
        viewer.plot()
        labelViewer.plot()


# In[ ]:





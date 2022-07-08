import sys

def package_info():
    packages = {}

    packages['python'] = sys.version.replace('\n', '| ')

    for pkg in ['fipy', 'numpy', 'pysparse', 'scipy', 'matplotlib', 'mpi4py', 'petsc4py', 'pyamgx']:
        try:
            mod = __import__(pkg)

            packages[pkg] = mod.__version__
        except ImportError as e:
            packages[pkg] = 'not installed'

        except Exception as e:
            packages[pkg] = 'version check failed: {}'.format(e)

    ## PyTrilinos
    try:
        import PyTrilinos
        packages['PyTrilinos'] = PyTrilinos.version()
    except ImportError as e:
        packages['PyTrilinos'] = 'not installed'
    except Exception as e:
        packages['PyTrilinos'] = 'version check failed: {}'.format(e)

    ## Mayavi uses a non-standard approach for storing its version number.
    try:
        from mayavi.__version__ import __version__ as mayaviversion
        packages['mayavi'] = mayaviversion
    except ImportError as e:
        try:
            from enthought.mayavi.__version__ import __version__ as mayaviversion
            packages['mayavi'] = mayaviversion
        except ImportError as e:
            packages['mayavi'] = 'not installed'
    except Exception as e:
        packages['mayavi'] = 'version check failed: {}'.format(e)

    ## Gmsh version
    try:
        from fipy.meshes.gmshMesh import gmshVersion
        gmshversion = gmshVersion()
        if gmshversion is None:
            packages['gmsh'] = 'not installed'
        else:
            packages['gmsh'] = gmshversion
    except Exception as e:
        packages['gmsh'] = 'version check failed: {}'.format(e)

    try:
        from fipy.solvers import solver_suite
        packages['solver'] = solver_suite
    except Exception as e:
        packages['solver'] = str(e)

    return packages

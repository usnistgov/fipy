import sys

def package_info():
    packages = ['python']
    versions = [sys.version.replace('\n', '| ')]

    for pkg in ['fipy', 'numpy', 'pysparse', 'scipy', 'matplotlib', 'mpi4py', 'petsc4py', 'pyamgx']:
        packages.append(pkg)

        try:
            mod = __import__(pkg)

            versions.append(mod.__version__)
        except ImportError as e:
            versions.append('not installed')

        except Exception as e:
            versions.append('version check failed: {}'.format(e))

    ## PyTrilinos
    packages.append('PyTrilinos')
    try:
        import PyTrilinos
        versions.append(PyTrilinos.version())
    except ImportError as e:
        versions.append('not installed')
    except Exception as e:
        versions.append('version check failed: {}'.format(e))

    ## Mayavi uses a non-standard approach for storing its version number.
    packages.append('mayavi')
    try:
        from mayavi.__version__ import __version__ as mayaviversion
        versions.append(mayaviversion)
    except ImportError as e:
        try:
            from enthought.mayavi.__version__ import __version__ as mayaviversion
            versions.append(mayaviversion)
        except ImportError as e:
            versions.append('not installed')
    except Exception as e:
        versions.append('version check failed: {}'.format(e))

    ## Gmsh version
    packages.append('gmsh')
    try:
        from fipy.meshes.gmshMesh import gmshVersion
        gmshversion = gmshVersion()
        if gmshversion is None:
            versions.append('not installed')
        else:
            versions.append(gmshversion)
    except Exception as e:
        versions.append('version check failed: {}'.format(e))

    packages.append("solver")
    try:
        from fipy.solvers import solver
        versions.append(solver)
    except Exception as e:
        versions.append(str(e))

    return (packages, versions)

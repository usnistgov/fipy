import json
import platform
import subprocess
import sys

__all__ = ["conda_info", "pip_info", "package_info", "platform_info", "nix_info"]

def conda_info(conda="conda"):
    """Collect information about conda environment.

    Parameters
    ----------
    conda : str
        Name of conda executable (default: "conda").

    Returns
    -------
    dict
        Result of `conda info` and `conda env export` for active conda
        environment.
    """
    info = {}
    p = subprocess.Popen([conda, "info", "--json"], stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    stdout = stdout.decode('ascii')

    info["conda_info"] = json.loads(stdout)

    p = subprocess.Popen([conda, "env", "export",
                          "--prefix", info["conda_info"]["active_prefix"],
                          "--json"],
                         stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    stdout = stdout.decode('ascii')
    info["conda_env"] = json.loads(stdout)

    return info

def pip_info(python="python"):
    """Collect information about pip environment.

    Parameters
    ----------
    python : str
        Name of Python executable (default: "python").

    Returns
    -------
    list of dict
        Result of `pip list --format json`.
    """
    info = {}
    p = subprocess.Popen([python, "-m", "pip", "list", "--format", "json"], stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    stdout = stdout.decode('ascii')

    info["pip"] = json.loads(stdout)

    return info

def nix_info():
    """Collect information about nix environment.

    Returns
    -------
    dict
        Result of `nix derivation show .#fipy `.
    """
    info = {}

    # Better output would be from
    #
    #     $ nix-store -q --tree $(nix-store --realize $(nix eval --raw .#fipy.drvPath))
    #
    # However, this is difficult to execute and can require another build.
    # Also doesn't return json.

    p = subprocess.Popen(
        [
         "nix",
         "derivation",
         "show",
         ".#fipy"
        ], stdout=subprocess.PIPE)

    stdout, _ = p.communicate()
    stdout = stdout.decode('ascii')

    info["nix"] = json.loads(stdout)

    return info

def package_info():
    """Collect information about installed packages FiPy uses.

    Returns
    -------
    dict
        Versions of important Python packages.
    """
    packages = {}

    packages['python'] = sys.version #.replace('\n', '| ')

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

def platform_info():
    """Collect information about platform Python is running in.

    Returns
    -------
    dict
        Data extracted from `platform` package.
    """
    return {
        "architecture": platform.architecture(),
        "machine": platform.machine(),
        "node": platform.node(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "release": platform.release(),
        "system": platform.system(),
        "version": platform.version()
    }

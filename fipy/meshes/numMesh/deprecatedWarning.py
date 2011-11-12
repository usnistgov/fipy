
import warnings

def numMeshDeprecated():
    tmp = "Meshes in `fipy.meshes.numMesh` are deprecated. Use meshes from `fipy.meshes`."
    warnings.warn(tmp, DeprecationWarning, stacklevel=3)


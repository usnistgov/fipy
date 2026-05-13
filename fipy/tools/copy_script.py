from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import sys
import typer
from typing_extensions import Annotated

__all__ = ["app"]

# from https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
def import_from_path(module_name, file_path):
    spec = spec_from_file_location(module_name, file_path)
    module = module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

def copy_script(
    copy_from: Annotated[str, typer.Argument(help="path and file name containing script to copy")],
    copy_to: Annotated[str, typer.Argument(help="path and file name to save script to")]
):
    copy_from = Path(copy_from)
    copy_to = Path(copy_to)

    if copy_to.exists():
        raise FileExistsError(copy_to)

    import fipy.tests.doctestPlus

    mod = import_from_path("copy_script_module", copy_from)
    script = fipy.tests.doctestPlus._getScript(name="copy_script_module")
    script = f"""
## This script was derived from
## '{copy_from}'

{script}
"""
    with open(copy_to, "w") as f:
        f.write(script)

    print(f"Script code exported from '{copy_from}' to '{copy_to}'")

app = typer.Typer()
app.command()(copy_script)

if __name__ == "__main__":
    app()

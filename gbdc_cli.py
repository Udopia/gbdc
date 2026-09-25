"""Console-script launcher for the bundled ``gbdc`` command-line tool.

The compiled CLI binary is shipped as package data next to this module (under
``gbdc.scripts/``) instead of in the wheel's scripts directory. Placing an ELF
executable in the scripts directory makes auditwheel relocate it and leave a
shim that resolves the binary through ``sysconfig.get_path("platlib")``, which
points at the global site-packages and therefore breaks ``pip install --user``
installs (pypa/auditwheel#340). This thin launcher instead locates the binary
relative to its own install location, so it works for every install layout.
"""

import os
import sys


def _binary_path() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    name = "gbdc.exe" if os.name == "nt" else "gbdc"
    return os.path.join(here, "gbdc.scripts", name)


def main() -> None:
    exe = _binary_path()
    args = [exe, *sys.argv[1:]]
    if os.name == "nt":
        import subprocess

        sys.exit(subprocess.call(args))
    os.execv(exe, args)


if __name__ == "__main__":
    main()

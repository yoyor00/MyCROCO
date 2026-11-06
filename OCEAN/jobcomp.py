#!/usr/bin/env python3
"""CROCO compilation driver (``jobcomp.py``, symlinked as ``jobcomp``).

    ./jobcomp [OPTIONS]              build CROCO in ./Compile
    ./jobcomp --sh [OPTIONS]         run the historical bash script jobcomp.sh
    ./jobcomp --help                 all the options

The driver itself lives in ``python/jobcomp`` (see its docstring); this file
only puts ``python/`` on ``sys.path`` and calls
``jobcomp.jobcomp_driver.main``.  ``jobcomp.conf`` and ``jobcomp.sh`` are
looked for in this directory.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "python"))

try:
    from jobcomp.jobcomp_driver import main
except ImportError:  # pragma: no cover - a broken/incomplete installation
    sys.stderr.write(
        f"ERROR: cannot import the 'jobcomp' package from {HERE / 'python'}\n"
    )
    raise SystemExit(1)

if __name__ == "__main__":
    sys.exit(main(script_path=__file__))

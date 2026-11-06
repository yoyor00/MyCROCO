"""Configuration: defaults, the optional ``jobcomp.conf`` and the lookup.

Configuration is resolved with, in decreasing priority:

1. command line options                  (``--src``, ``--fc``, ...)
2. ``CROCO_*`` environment variables     (``CROCO_SRC``, ``CROCO_CFT1``, ...)
3. a ``jobcomp.conf`` file next to ``jobcomp.py``
4. the defaults of :class:`Config` (the "user's modifications" block of the
   bash script)

``jobcomp.conf`` holds one ``KEY=VALUE`` per line using the historical names,
which is how ``create_config.bash`` configures a generated ``CROCO_IN``: no
more ``sed``-ing of the script itself.  Values are expanded (``$VAR``, ``~``)
when they are used, so a file may contain ``SOURCE1=${OCE}`` and stay valid on
any machine.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Optional

from .errors import BuildError


class Config:
    """The settings of one run: where the script is, and what it should use."""

    CONF_NAME = "jobcomp.conf"
    LEGACY_SCRIPT = "jobcomp.sh"  # the historical bash implementation

    #: ``CROCO_*`` environment variable -> conf key / default attribute.
    CONF_KEYS = {
        "CROCO_SRC": "SOURCE1",
        "CROCO_CFT1": "FC",
        "CROCO_MPIF90": "MPIF90",
        "CROCO_MPILIB": "MPILIB",
        "CROCO_MPIINC": "MPIINC",
        "CROCO_PRISM_ROOT_DIR": "PRISM_ROOT_DIR",
        "CROCO_XIOS_ROOT_DIR": "XIOS_ROOT_DIR",
    }

    # -----------------------------------------------------------------------
    # BEGIN OF USER'S MODIFICATIONS
    #
    # Defaults only: prefer a jobcomp.conf file or the environment so that the
    # driver itself never has to be edited.
    # -----------------------------------------------------------------------
    SOURCE1 = "../croco/OCEAN"
    FC = "gfortran"
    MPIF90 = "mpif90"
    MPILIB = ""
    MPIINC = ""
    PRISM_ROOT_DIR = "../../../oasis3-mct/compile_oa3-mct"
    XIOS_ROOT_DIR = "$HOME/xios"
    # -----------------------------------------------------------------------
    # END OF USER'S MODIFICATIONS
    # -----------------------------------------------------------------------

    def __init__(
        self,
        script_path: Optional[str] = None,
        environ: Optional[Dict[str, str]] = None,
    ) -> None:
        #: Real path of ``jobcomp.py``: a checkout and an installed (copied)
        #: driver then both find their ``jobcomp.conf`` and the bash
        #: ``jobcomp`` next to the real file.
        self.script_path = (
            os.path.realpath(str(script_path)) if script_path else None
        )
        self.environ = os.environ if environ is None else environ

    # -- files next to the driver -----------------------------------------
    @property
    def script_dir(self) -> Path:
        """Directory holding ``jobcomp.py``, ``jobcomp.sh`` and its conf."""
        if self.script_path is not None:
            return Path(self.script_path).parent
        return Path(os.path.realpath(sys.argv[0])).parent

    def load_conf(self, path: Optional[Path] = None) -> Dict[str, str]:
        """Read the optional ``jobcomp.conf`` (next to ``jobcomp.py``)."""
        try:
            text = Path(path or self.script_dir / self.CONF_NAME).read_text()
        except OSError:
            return {}
        values: Dict[str, str] = {}
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            value = value.strip()
            if len(value) > 1 and value[0] == value[-1] and value[0] in "\"'":
                value = value[1:-1]
            values[key.strip()] = value
        return values

    # -- values ------------------------------------------------------------
    @staticmethod
    def expand(path: str) -> Path:
        """Turn a configured path into a ``Path`` (``~``/``$VAR`` resolved)."""
        return Path(os.path.expanduser(os.path.expandvars(str(path))))

    def default(self, env_name: str) -> str:
        """The hard-coded default behind ``env_name`` (``""`` if unknown)."""
        return getattr(self, self.CONF_KEYS.get(env_name, ""), "")

    def pick(
        self, cli_value, env_name: str, default: Optional[str] = None
    ) -> str:
        """CLI > environment > ``jobcomp.conf`` > hard-coded default."""
        if cli_value is not None:
            return str(cli_value)
        value = self.environ.get(env_name)
        if value:
            return value
        conf_value = self.load_conf().get(self.CONF_KEYS.get(env_name, ""), "")
        if conf_value:
            return os.path.expanduser(os.path.expandvars(conf_value))
        return self.default(env_name) if default is None else default

    def legacy_script(self) -> Path:
        """Locate ``jobcomp``, the unmodified bash implementation."""
        script = self.script_dir / self.LEGACY_SCRIPT
        if not script.is_file():
            raise BuildError(f"cannot find {self.LEGACY_SCRIPT}")
        return script

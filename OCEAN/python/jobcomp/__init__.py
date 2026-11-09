"""CROCO compilation driver: the Python successor of ``OCEAN/jobcomp``.

    ./jobcomp [OPTIONS]              build CROCO in ./Compile
    ./jobcomp --sh [OPTIONS]         run the historical bash script jobcomp.sh
    ./jobcomp --help                 all the options

Settings, highest priority first: command line, ``CROCO_*`` environment
variables, the ``jobcomp.conf`` file next to it, and the defaults
in :class:`jobcomp.config.Config`.  The produced ``Compile/`` tree,
``Makedefs`` and ``CROCOFILES`` headers are the same as the bash script's.

Another build system (scons, cmake, ...) can read the resolved flags and
paths with ``Jobcomp(...).parameters()`` -- see :mod:`jobcomp.parameters`.

Modules (all pure libraries; ``OCEAN/jobcomp.py`` is the entry point):

    errors, config, settings, shell, copy, toolchain, preprocessor, features,
    context, makedefs, build, parameters, jobcomp_driver -- the program,
    one class per module
    components/           -- the build features, one module per feature

Requires Python 3.8+ and the standard library only.
"""

__version__ = "1.0"

"""CROCO compilation driver: the Python successor of ``OCEAN/jobcomp``.

    ./jobcomp [OPTIONS]              build CROCO in ./Compile
    ./jobcomp --scons [OPTIONS]      build with the optional SCons backend
    ./jobcomp --sh [OPTIONS]         run the historical bash script jobcomp.sh
    ./jobcomp --help                 all the options

Settings, highest priority first: command line, ``CROCO_*`` environment
variables, the ``jobcomp.conf`` file next to it, and the defaults
in :class:`jobcomp.config.Config`.  The produced ``Compile/`` tree,
``Makedefs`` and ``CROCOFILES`` headers are the same as the bash script's.
The SCons backend also reads ``Compile/croco_build.json``, written from the
resolved build context just before compilation.

Another build system (scons, cmake, ...) can read the resolved flags and
paths with ``Jobcomp(...).parameters()`` -- see :mod:`jobcomp.parameters`.

Modules (all pure libraries; ``OCEAN/jobcomp.py`` is the entry point):

    errors, config, settings, shell, copy, flags, makefile, manifest,
    toolchain, preprocessor, features, context, makedefs, build, parameters,
    jobcomp_driver -- the program, one class per module
    components/           -- the build features, one module per feature
    ../scons/             -- the optional SCons compilation backend

Requires Python 3.8+ and the standard library only.  The optional SCons
backend additionally requires the third-party ``scons`` package.
"""

__version__ = "1.0"

#!/usr/bin/env python3
##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""
This script aimed at benchmarking CROCO to compare the various parallel modes.
"""

##########################################################
import os
import sys

# This script only writes figures to files (PNG/SVG), it never opens an
# interactive window. On some clusters DISPLAY is set (X forwarding) but the
# X server does not support XInput2, which makes the GTK backend abort with
# "Gdk-ERROR **: XInput2 support not found on display" as soon as pyplot is
# imported. Force the non-interactive Agg backend (still overridable by the
# user through MPLBACKEND).
os.environ.setdefault("MPLBACKEND", "Agg")

from libs.config import Config
from libs.helpers import print_exception
from libs.benchmarking import Benchmarking


##########################################################
def main():
    # build config
    config = Config()
    config.parse()
    config.load_config()

    if config.globalhtml:
        config.html_global_report()
    else:
        # benchmark
        bench = Benchmarking(config)
        bench.run()


##########################################################
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print_exception(e)
        sys.exit(1)

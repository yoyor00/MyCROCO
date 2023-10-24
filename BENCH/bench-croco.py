#!/usr/bin/env python3
##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
This script aimed at benchmarking CROCO to compare the various parallel modes.
'''

##########################################################
from libs.config import Config
from libs.benchmarking import Benchmarking

##########################################################
def main():
    # build config
    config = Config()
    config.parse()
    config.load_config()

    # benchmark
    bench = Benchmarking(config)
    bench.run()

##########################################################
if __name__ == '__main__':
    main()

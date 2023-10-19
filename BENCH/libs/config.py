##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json5
import argparse

##########################################################
class Config:
    '''Configure the benchmarking'''
    def __init__(self):
        self.args = None

    def parse(self) -> None:
        '''Parse the program options to init the run'''
        parser = argparse.ArgumentParser(
            prog = 'bench-croco.py',
            description = 'Automate the benchmarking of CROCO by doing rewrite.')

        # register options
        parser.add_argument('-w', '--workdir', help="Working sub-directory name.", default='rundir')
        parser.add_argument('-c', '--cases', help="Specific case to run.", default='@default')
        parser.add_argument('-v', '--variants', help="Variant to run.", default='@default')
        parser.add_argument('-V', '--verbose', help="Enable the full logging instead of printing just summarys.", action='store_true')
        parser.add_argument('-r', '--rebuild', help="Force rebuilding.", action='store_true')
        parser.add_argument('-s', '--steps', help="Number of time steps to perform.")
        parser.add_argument('--check', help="Check that the resulting mesh is same than seq run.", action='store_true')
        parser.add_argument('-C', '--config', help="Config file to use.", default='config.jsonc')
        parser.add_argument('-m', '--modes', help="List of modes to run.", default='build,run,plot')
        parser.add_argument('-R', '--runs', help="Number or runs to perform.", default='4')
        parser.add_argument('-j', '--jobs', help="Make -j option value to build", default='8')
        parser.add_argument('--results', help="Name of the results directory to use.", default='results')

        # parse
        self.args = parser.parse_args()

        # reshape
        self.modes = self.args.modes.split(',')
        self.variant_names = self.args.variants.split(',')
        self.case_names = self.args.cases.split(',')
        self.verbose = self.args.verbose
        self.workdir = self.args.workdir
        self.steps = self.args.steps
        self.rebuild = self.args.rebuild
        self.runs = int(self.args.runs)
        self.make_jobs = int(self.args.jobs)
        self.results = self.args.results

        # extract some many time used paths
        self.croco_source_dir = os.path.abspath(f"{__file__}/../../../")

    def load_config(self) -> None:
        '''Load the config file and extract some infos'''
        # load
        with open(self.args.config, 'r') as fp:
            self.config = json5.load(fp)

        # apply default
        if self.variant_names == ["@default"]:
            self.variant_names = self.config['default_selection']['variants']
        if self.case_names == ['@default']:
            self.case_names = self.config['default_selection']['cases']

        # apply user selection
        if self.variant_names == ['@all']:
            self.variant_names = self.config['variants'].keys()
        if self.case_names == ['@all']:
            self.case_names = self.config['cases'].keys()

        # extract variant groups
        final_variants = []
        for variant in self.variant_names:
            if variant.startswith('@'):
                final_variants += self.config['meta_variants'][variant]
            else:
                final_variants.append(variant)
        self.variant_names = final_variants

##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json5
import yaml
import glob
from yaml import Loader
from jsonschema import validate
import argparse
import platform
import copy
from datetime import datetime
from .messaging import Messaging
from .configfile import ConfigFile

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
        parser.add_argument('-t', '--title', help="A possible extra name to prepent to the result directory name", default='.')
        parser.add_argument('-a', '--auto-skip', help='Automatically skip the cases with not enougth ressources.', action='store_true')
        parser.add_argument('-N', '--no-previous', help='Do not load the preivous missing results to plot full graph.', action='store_true')

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
        self.title = self.args.title
        self.auto_skip = self.args.auto_skip
        self.no_previous = self.args.no_previous

        # compute clean result subdir name
        hostname = platform.node()
        run_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.results = os.path.join(self.args.results, self.title, f"{hostname}-{run_date}")

        # pattern to search results files (to also take the previous runs if not re-run all)
        if self.no_previous:
            self.results = self.results
        else:
            self.results_pattern = os.path.join(self.args.results, self.title, f"{hostname}-*")

        # extract some many time used paths
        self.croco_source_dir = os.path.abspath(f"{__file__}/../../../")

    def filter_variant_ressources(self):
        # extract needed vars
        variants = self.config['variants']
        host = self.host
        host_ressources = host['ressources']

        # filter
        for name in self.variant_names.copy():
            variant = variants[name]
            for ressource_name, ressource_min in variant['requires'].items():
                has = host_ressources[ressource_name]
                if has < int(ressource_min):
                    Messaging.step(f"Filter out {name} due to lack of {ressource_name} ({has} < {ressource_min})...")
                    if self.auto_skip:
                        self.variant_names.remove(name)
                    else:
                        raise Exception(f"Cannot run {name} due to lack of {ressource_name} ({has} < {ressource_min})... You can use -a/--auto-skip to skip silently.")

    def load_config(self) -> None:
        '''Load the config file and extract some infos'''
        # message
        Messaging.section("Loading configuration")

        # load
        self.config = ConfigFile(self.args.config).load()

        # apply meta : @.....
        self.apply_meta()

        # select host
        self.select_host()

        # filter ressources
        self.filter_variant_ressources()

    def apply_meta(self):
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

    def select_host(self):
        # select tunning
        hostname = platform.node()
        hosts = self.config['hosts']
        if hostname in hosts:
            self.host = hosts[hostname]
        else:
            self.host = hosts['@default']

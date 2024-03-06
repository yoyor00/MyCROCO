##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import argparse
import platform
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

        # get hostname from system
        selected_host = platform.node()

        # register options
        parser.add_argument('-w', '--workdir', help="Working sub-directory name.", default='rundir')
        parser.add_argument('-c', '--cases', help="Specific case to run.", default='@default')
        parser.add_argument('-v', '--variants', help="Variant to run.", default='@default')
        parser.add_argument('-V', '--verbose', help="Enable the full logging instead of printing just summaries.", action='store_true')
        parser.add_argument('-r', '--rebuild', help="Force rebuilding.", action='store_true')
        parser.add_argument('--check', help="Check that the resulting mesh is same than seq run.", action='store_true')
        parser.add_argument('-C', '--config', help="Config file to use.", default='config.jsonc')
        parser.add_argument('-m', '--modes', help="List of modes to run.", default='build,run,check,plot')
        parser.add_argument('-R', '--runs', help="Number or runs to perform.", default='4')
        parser.add_argument('-j', '--jobs', help="Make -j option value to build", default='8')
        parser.add_argument('--results', help="Name of the results directory to use.", default='results')
        parser.add_argument('-t', '--title', help="A possible extra name to prepent to the result directory name", default='.')
        parser.add_argument('-a', '--auto-skip', help='Automatically skip the cases with not enougth ressources.', action='store_true')
        parser.add_argument('-N', '--no-previous', help='Do not load the preivous missing results to plot full graph.', action='store_true')
        parser.add_argument(      '--build-ref', help="Build the reference directory to be stored once somewhere for validation.", type=str, default=False)
        parser.add_argument('-u', '--use-ref', help="Use the reference directory as source of compare instead of seq run.", type=str, default=False)
        parser.add_argument(      '--jobcomp', help="Force using jobcomp instead of cmake.", action='store_true')
        parser.add_argument(      '--host', help="Force host config to use.", type=str, default=selected_host)
        parser.add_argument(      '--enable-debug', help="Use cmake debug build instead of release.", action='store_true')
        parser.add_argument(      '--rvtk', help="Enable usage of RVTK_DEBUG and give name of variant reference.", action='store_true')
        parser.add_argument(      '--compare-to', help="Set which variant is used as reference (default='sequential').", type=str, default='sequential')

        # parse
        self.args = parser.parse_args()

        # reshape
        self.modes = self.args.modes.split(',')
        self.variant_names = self.args.variants.split(',')
        self.case_names = self.args.cases.split(',')
        self.verbose = self.args.verbose
        self.capture = not self.verbose
        self.workdir = os.path.abspath(self.args.workdir)
        self.rebuild = self.args.rebuild
        self.runs = int(self.args.runs)
        self.make_jobs = int(self.args.jobs)
        self.title = self.args.title
        self.auto_skip = self.args.auto_skip
        self.no_previous = self.args.no_previous
        self.build_ref = self.args.build_ref
        self.use_ref = self.args.use_ref
        self.force_jobcomp = self.args.jobcomp
        self.use_host_config = self.args.host
        self.debug_build = self.args.enable_debug
        self.variant_ref_name = self.args.compare_to
        self.rvtk = self.args.rvtk

        # compute clean result subdir name
        use_host_config = self.use_host_config
        run_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.results = os.path.join(self.args.results, self.title, f"{use_host_config}-{run_date}")

        # pattern to search results files (to also take the previous runs if not re-run all)
        if self.no_previous:
            self.results = self.results
        else:
            self.results_pattern = os.path.join(self.args.results, self.title, f"{use_host_config}-*")

        # extract some many time used paths
        self.croco_source_dir = os.path.abspath(f"{__file__}/../../../")

        # detect old way or cmake way
        self.has_cmake = os.path.exists(os.path.join(self.croco_source_dir, 'CMakeLists.txt'))

        # detect if is minicroco jobcomp
        self.is_minicroco_jobcomp = os.path.exists(os.path.join(self.croco_source_dir, 'OCEAN/config/cppdefs.h'))

    def pre_checks(self):
        # basic check
        if self.use_ref and not os.path.exists(self.use_ref):
            raise Exception(f"You gave --use-ref={self.use_ref}, but directory does not exist !")
        if self.build_ref and not self.variant_ref_name in self.variant_names:
            raise Exception(f"You gave --build-ref={self.use_ref}, but not using the required '{self.variant_ref_name}' variant !")

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
        # load
        self.config = ConfigFile(self.args.config).load()
        
        # message
        Messaging.section("Applying selection")

        # apply meta : @.....
        self.apply_meta()

        # select host
        self.select_host()

        # filter ressources
        self.filter_variant_ressources()

        # check
        self.pre_checks()

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

        # extract cases groups
        final_cases = []
        for case in self.case_names:
            if case.startswith('@'):
                final_cases += self.config['meta_cases'][variant]
            else:
                final_cases.append(case)
        self.case_names = final_cases

    def select_host(self):
        # get it
        hostname = platform.node()

        # select tunning
        selected_host = self.use_host_config
        hosts = self.config['hosts']
        
        # fallback
        if not selected_host in hosts:
            selected_host = '@default'
            
        # display
        Messaging.step(f"Hostname : {hostname}")
        Messaging.step(f"Select host : {selected_host}")
        
        # apply
        self.host = hosts[selected_host]

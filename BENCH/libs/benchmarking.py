##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json
import subprocess
from typing import Union
from .config import Config
from .croco import Croco
from .messaging import Messaging
from .system import gen_system_info
from .plotting import Plotting
from .helpers import run_shell_command

##########################################################
class Benchmarking:
    def __init__(self, config: Config):
        # set
        self.config = config

        # gen str
        result_dir = self.config.results
        case_names = ', '.join(self.config.case_names)
        variant_names = ', '.join(self.config.variant_names)

        # init
        Messaging.section("Benchmark initialization")
        Messaging.step(f"Results  = {result_dir}")
        Messaging.step(f"Cases    = {case_names}")
        Messaging.step(f"Variants = {variant_names}")

        # create work dir
        os.makedirs(config.workdir, exist_ok=True)

        # create result dir
        os.makedirs(config.results, exist_ok=True)

        # create last directory
        self.create_last_symlink()

        # make sequential first as we use as reference to validate the data produced
        self.make_ref_variant_first()

        # Create croco instances
        self.instances = self.create_croco_instances()

        # dump
        self.dump_bench_infos()

    def create_last_symlink(self):
        # extract
        results = self.config.results

        # create last directory
        dirname = os.path.dirname(results)
        basename = os.path.basename(results)
        link_name = f"{dirname}/last"
        if os.path.exists(link_name):
            os.remove(link_name)
        os.symlink(basename, link_name)

    def make_ref_variant_first(self):
        '''
        Make sequential first as we use as reference to validate the data produced
        '''
        # get ref
        ref_name = self.config.variant_ref_name

        # ref variant
        if ref_name in self.config.case_names:
            self.config.variant_names.remove(ref_name)
            self.config.variant_names.insert(0, ref_name)

    def create_croco_instances(self) -> Union[Croco]:
        # extract some
        config = self.config

        # res
        res = []

        # build combination
        for case_name in config.case_names:
            for variant_name in config.variant_names:
                # get config to know how to filter
                case_config = config.config['cases'][case_name]

                # filter
                if not variant_name in case_config.get('unsupported', []):
                    res.append(Croco(config, case_name, variant_name))
                else:
                    Messaging.step(f"Skip unsupported {case_name}/{variant_name}")

        # ok
        return res

    def run(self):
        # vars
        cnt = len(self.instances)
        variant_ref_name = self.config.variant_ref_name

        # build
        if 'build' in self.config.modes:
            for id, instance in enumerate(self.instances):
                instance.build(extra_info=f" - [ {id + 1} / {cnt} ]", force_rebuild = self.config.rebuild)

        # run
        if 'run' in self.config.modes:
            instance: Croco
            for id, instance in enumerate(self.instances):
                # run
                instance.run(extra_info=f" - [ {id + 1} / {cnt} ]")
                # check
                if 'check' in self.config.modes:
                    instance.check()
                # if need to store the ref
                if self.config.build_ref and instance.variant_name == variant_ref_name:
                    instance.make_ref()

        # plot
        if 'plot' in self.config.modes:
            self.plot()

        # report
        if self.config.report:
            self.config.report.save()
            self.config.report.print()

    def dump_bench_infos(self):
        # dump bench infos
        Messaging.section("Dumping system infos")

        # get needs
        results = self.config.results
        
        # gen
        system = gen_system_info()
        with open(f"{results}/system.json", "w+") as fp:
            json.dump(system, fp=fp, indent='\t')

        # dump config
        with open(f"{results}/config.json", "w+") as fp:
            json.dump(self.config.config, fp=fp, indent='\t')

        # display some
        hostname = system['plateform']['hostname']
        processor_name = system['plateform']['processor_name']
        Messaging.step(f"Hostname  : {hostname}")
        Messaging.step(f"Processor : {processor_name}")

        # dump the CPU infos
        run_shell_command(f"hwloc-ls --of console {results}/cpu.txt")
        run_shell_command(f"hwloc-ls --of svg {results}/cpu.svg")

    def plot(self):
        plot = Plotting(self.config)
        plot.plot()

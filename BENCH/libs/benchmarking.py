##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json
from .config import Config
from .croco import Croco
from .messaging import Messaging
from .system import gen_system_info
from .plotting import Plotting

##########################################################
class Benchmarking:
    def __init__(self, config: Config):
        # set
        self.config = config

        # init
        Messaging.section("Benchmark initialization")
        Messaging.step(f"Results  = {self.config.results}")
        Messaging.step(f"Cases    = {self.config.case_names}")
        Messaging.step(f"Variants = {self.config.variant_names}")

        # create work dir
        os.makedirs(config.workdir, exist_ok=True)

        # create result dir
        os.makedirs(config.results, exist_ok=True)

        # create last directory
        self.create_last_symlink()

        # make sequential first as we use as reference to validate the data produced
        self.make_sequential_first()

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

    def make_sequential_first(self):
        '''
        Make sequential first as we use as reference to validate the data produced
        '''
        if 'sequential' in self.config.case_names:
            self.config.variant_names.remove('sequential')
            self.config.variant_names.insert(0, 'sequential')

    def create_croco_instances(self) -> [Croco]:
        # res
        res = []

        # build combination
        for case in self.config.case_names:
            for variant in self.config.variant_names:
                res.append(Croco(self.config, case, variant))

        # ok
        return res

    def run(self):
        cnt = len(self.instances)
        # build
        if 'build' in self.config.modes:
            for id, instance in enumerate(self.instances):
                instance.build(extra_info=f" - [ {id + 1} / {cnt} ]", force_rebuild = self.config.rebuild)
        # run
        if 'run' in self.config.modes:
            for id, instance in enumerate(self.instances):
                # run
                instance.run(extra_info=f" - [ {id + 1} / {cnt} ]")
                # check
                if 'check' in self.config.modes:
                    instance.check()

        # plot
        if 'plot' in self.config.modes:
            self.plot()

    def dump_bench_infos(self):
        # get needs
        results = self.config.results
        
        # gen
        system = gen_system_info()
        with open(f"{results}/system.json", "w+") as fp:
            json.dump(system, fp=fp, indent='\t')

        # dump config
        with open(f"{results}/config.json", "w+") as fp:
            json.dump(self.config.config, fp=fp, indent='\t')

    def plot(self):
        plot = Plotting(self.config)
        plot.plot()

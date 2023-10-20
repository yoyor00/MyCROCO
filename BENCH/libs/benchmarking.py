##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
from .config import Config
from .croco import Croco
from .messaging import Messaging

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

        # make sequential first as we use as reference to validate the data produced
        self.make_sequential_first()

        # Create croco instances
        self.instances = self.create_croco_instances()

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
        for id, instance in enumerate(self.instances):
            instance.build(extra_info=f" - [ {id + 1} / {cnt} ]")
        # run
        for id, instance in enumerate(self.instances):
            instance.run(extra_info=f" - [ {id + 1} / {cnt} ]")

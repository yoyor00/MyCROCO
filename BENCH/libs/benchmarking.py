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

        # create work dir
        if not os.path.exists(config.workdir):
            os.mkdir(config.workdir)

        # create result dir
        if not os.path.exists(config.results):
            os.mkdir(config.results)

        # make sequential first as we use as reference to validate the data produced
        self.make_sequential_first()

        # Create croco instances
        self.instances = self.create_croco_instances()

    def make_sequential_first(self):
        '''
        Make sequential first as we use as reference to validate the data produced
        '''
        if 'sequential' in self.config.case_names:
            self.config.case_names.removes('sequential')
            self.config.case_names.insert(0, 'sequential')

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
        for instance in self.instances:
            instance.build()

##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import glob
import json
import numpy
from matplotlib import pyplot
from .config import Config
from .messaging import Messaging

##########################################################
class Plotting:
    def __init__(self, config: Config):
        self.config = config

    def load_case_variant_data(self, case_name: str, variant_name: str) -> dict:
        # name
        name = f"{case_name}-{variant_name}"
        results_pattern = self.config.results_pattern

        # calc
        fname = f"result-{name}.json"
        avail = glob.glob(f"{results_pattern}/{fname}")

        # has none
        if avail is None or len(avail) == 0:
            return None
        
        # we can sort by date and take the last one
        avail.sort()
        fpath = avail[-1]

        # progress
        Messaging.step(f"Loading {fpath}")

        # load it
        with open(fpath, 'r') as fp:
            return json.load(fp)

    def load_variants(self, data: dict, case_name: str) -> bool:
        # extract vars
        variant_names = self.config.config['variants'].keys()

        # loop on all variants
        for variant_name in variant_names:
            # load json
            case_variant_data = self.load_case_variant_data(case_name, variant_name)

            # merge
            if case_variant_data is not None:
                data['variants'].append(variant_name)
                data['means'].append(case_variant_data['results'][0]['mean'])
                data['median'].append(case_variant_data['results'][0]['median'])
                data['stddev'].append(case_variant_data['results'][0]['stddev'])
                data['min'].append(case_variant_data['results'][0]['min'])
                data['max'].append(case_variant_data['results'][0]['max'])
                data['q1'].append(numpy.quantile(case_variant_data['results'][0]['times'], 0.25))
                data['q3'].append(numpy.quantile(case_variant_data['results'][0]['times'], 0.75))

        # ok
        return True

    def build_data(self):
        # extract vars
        #case_names = self.config.case_names
        case_names = self.config.config['cases'].keys()

        # to fill
        data = {}

        # loop
        for name in case_names:
            # build
            case_data = {
                'variants': [],
                'means': [],
                'median': [],
                'min': [],
                'max': [],
                'q1': [],
                'q3': [],
                'stddev': [],
            }

            # load
            self.load_variants(case_data, name)

            # attach
            if len(case_data['variants']) != 0:
                data[name] = case_data

        # return
        return data

    def plot(self):
        # info
        Messaging.section(f'Plotting')


        # extract
        cases = self.config.config['cases']

        # get data
        data = self.build_data()

        # plot
        for case_name, case_data in data.items():
            # progress
            Messaging.step(f"Plotting {case_name}")

            # extract title
            case_config = cases[case_name]
            results = self.config.results
            title = case_config.get('title', case_name)
            

            # start plot
            fig, ax = pyplot.subplots()

            # fix some names
            # TODO : remove
            for key,value in enumerate(case_data['variants']):
                case_data['variants'][key] = value.replace('-wip',' progress')

            # build graph
            x_pos = numpy.arange(len(case_data['variants']))
            #ax.bar(x_pos, entry['median'], yerr=[entry['min'],entry['max']], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax.bar(x_pos, case_data['median'], yerr=case_data['stddev'], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax.set_ylabel('Average runtime (seconds)')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(case_data['variants'], rotation=45, ha='right')
            ax.set_title(f'CROCO - {title}')
            ax.yaxis.grid(True)

            # Save the figure and show
            pyplot.tight_layout()
            pyplot.savefig(f'{results}/plot-{case_name}.png')
            pyplot.savefig(f'{results}/plot-{case_name}.svg')
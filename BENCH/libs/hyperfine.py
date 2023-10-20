
##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Hyperfie (https://github.com/sharkdp/hyperfine) is a runner to measure perf
of a command.
I provide here a python short replacement as it is not installed by default
systems.
'''

##########################################################
import numpy
import json
import shutil
from tqdm import tqdm
from tempfile import NamedTemporaryFile
from .messaging import Messaging
from .helpes import run_shell_command_time, run_shell_command

##########################################################
def emulate_hyperfine(command: str, progress_bar: tqdm = None, runs: int = 2, retries: int = 8, verbose: bool = False) -> dict:
        # prepare some vars
        measures = []

        # display
        Messaging.command(command)

        # update progress bar
        if progress_bar is not None:
            progress_bar.refresh()

        # repeat
        for iteration in range(runs):
            # run
            Messaging.step(f"Running iteration - [ {iteration + 1} / {runs} ]")
            runtime = run_shell_command_time(f"{command}", verbose)
            Messaging.step(f"Runtime = {runtime}")
            measures.append(runtime)

            # progress forward
            if progress_bar is not None:
                progress_bar.update(1)

        # build summary
        Messaging.step(f"All = {measures}")
        summary = {
            'results': [
                {
                    'command': command,
                    'mean': numpy.mean(measures),
                    'median': numpy.median(measures),
                    'stddev': numpy.std(measures),
                    'min': numpy.min(measures),
                    'max': numpy.max(measures),
                    'times': measures
                }
            ]
        }

        # return
        return summary

##########################################################
def native_hyperfine(command: str, runs: int = 2) -> dict:
    with NamedTemporaryFile() as fp:
        hyper_command = f"hyperfine --export-json={fp.name} --warmup 0 --output /tmp/croco.log --runs {runs} \"{command}\""
        try:
            run_shell_command(hyper_command, capture=False)
            result = json.load(fp)
            return result
        except Exception as e:
            with open('/tmp/croco.log', 'r') as fp_log:
                content = fp_log.read()
                print(content)
            raise Exception("CROCO failed to run !")

##########################################################
def run_hyperfine(command: str, progress_bar: tqdm = None, runs: int = 2, retries: int = 8):
    hyperfine = False #shutil.which("hyperfine")
    if hyperfine:
        return native_hyperfine(command, runs)
    else:
        return emulate_hyperfine(command, progress_bar, runs)

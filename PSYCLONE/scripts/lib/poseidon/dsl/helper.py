##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import json
# psyclone
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.nodes import Call
# internal
from ..dsl.kernel import Kernel, KernelList, Loop, Routine
from ..base.psyir_helpers import get_first_loop, get_previous
from ..parsing.walker import Walker
from ..parsing.kernel_extractor import KernelExtractor

##########################################################
def get_psyir_from_code(code: str, free_form=False):
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(code, free_form=free_form)
    return psyir_tree

##########################################################
def extract_kernel_no_walk(code: str, free_form=False) -> Kernel:
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(code, free_form=free_form)
    loop = get_first_loop(psyir_tree)
    kernel = Kernel(loop)
    return kernel

##########################################################
def extract_kernels(code: str, free_form=False) -> KernelList:
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(code, free_form=free_form)
    return extract_kernels_from_psyir(psyir_tree)

##########################################################
def extract_kernels_from_file(code: str, free_form=False) -> KernelList:
    reader = FortranReader()
    psyir_tree = reader.psyir_from_file(code, free_form=free_form)
    return extract_kernels_from_psyir(psyir_tree)

##########################################################
def extract_kernels_from_psyir_internal(kernels, psyir_tree, ignore_loops=[]):
    top_loop: Loop
    for top_loop in psyir_tree.walk(Loop, stop_type=Loop):
        if top_loop.variable.name in ignore_loops:
            extract_kernels_from_psyir_internal(kernels, top_loop.loop_body, ignore_loops)
        else:
            kernels.push_kernel(Kernel(top_loop))

##########################################################
def extract_kernels_from_psyir(psyir_tree, ignore_loops=[]) -> KernelList:
    kernels = KernelList()
    extract_kernels_from_psyir_internal(kernels, psyir_tree, ignore_loops=ignore_loops)
    return kernels

##########################################################
def extract_kernel_walk(code: str, free_form=False) -> Kernel:
    return extract_kernels(code, free_form=free_form).get(0)

##########################################################
def dump_async_var_streams(psy_ir_tree):
    # loop on all calls
    for call in psy_ir_tree.walk(Call):
        previous_loop = get_previous(call, Loop)
        if previous_loop != None:
            fname = f"/tmp/{call.routine.name}.json"
            with open(fname, 'w+') as fp:
                json.dump(previous_loop.aync_var_last_stream, fp)

    # save state at the end of routing
    for routine in psy_ir_tree.walk(Routine):
        fname = f"/tmp/{routine.name}_out.json"
        with open(fname, 'w+') as fp:
            json.dump(routine.sync_var_last_stream, fp)

##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Some helper function to make some specific transform on the CROCO code.

Remark
------
There is a lot of things currently hard coded and which requires to be made
generic instead.
'''

##########################################################
# intenral
from .loops_helpers import is_loop_using_var
from .directives.ACCSetDeviceNumDirective import ACCSetDeviceNumDirective
# psyclone
from psyclone.psyir.nodes import Node, Routine, Literal
from psyclone.psyir.symbols import DataSymbol, INTEGER_TYPE, BOOLEAN_TYPE
from psyclone.psyir.nodes import Call, IntrinsicCall, Node, Loop, Reference
from psyclone.psyir.nodes.acc_directives import ACCLoopDirective
from psyclone.core import Signature
from psyclone.transformations import ACCLoopTrans, TransformationError, ACCEnterDataTrans
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir, InlinedKern
# poseidon
from poseidon.dsl.helper import KernelList, extract_kernels_from_psyir

##########################################################
def apply_acc_fetch_vars(psy) -> None:
    """
    Fetch vars used inside kernel loop.
    These vars can then be made private.

    TODO: rename function properly
    """
    vars_to_fetch = []
    for invoke in psy.invokes.invoke_list:
        for ref in invoke.schedule.walk(Reference):
            if ref.name.lower() in ['fc1d', 'cf1d', 'dc1d']:
                vars_to_fetch.append(ref.name.lower())
        vars_to_fetch = list(set(vars_to_fetch))

        if len(vars_to_fetch) > 0:
            # build sigs
            signatures = []
            for var in vars_to_fetch:
                signatures.append(Signature(var))

            ACCEnterDataTrans().apply(invoke.schedule, options={'signatures': signatures})
            for dir in invoke.schedule.walk(AccEnterDataDir):
                dir.sig_set = signatures

##########################################################
def apply_acc_loop_and_collapse(kernels: KernelList, options: dict) -> None:
    """
    For all kernel loops, collapse loop if possible.

    TODO: Make sure that this is always valid
    """
    for kernel in kernels.kernels:
        loop = kernel.root_node

        try:
            ACCLoopTrans().apply(loop, options=options)
        except TransformationError as err:
            # This loop can not be transformed, proceed to next loop
            print("Loop not parallelised because:", str(err))
            continue

        # Count the number of perfectly nested loops & collapse them
        num_nested_loops = 0
        next_loop = loop
        while isinstance(next_loop, Loop):
            num_nested_loops += 1
            if len(next_loop.loop_body.children) > 1:
                break
            next_loop = next_loop.loop_body.children[0]

        if num_nested_loops > 1:
            loop.parent.parent.collapse = num_nested_loops

##########################################################
def apply_acc_kernel(container: Node, collapse: bool, ignore_loops: list=[], merge_joinable: bool = False, options: dict = None) -> None:
    '''
    Apply the ACC kernel directives and make source optimization (fuse...)
    '''

    # extract kernels
    kernels = extract_kernels_from_psyir(container, ignore_loops=ignore_loops)

    # skip some kernels to tests
    # TODO fin a nice way to make this cleanly
    for kernel in kernels.kernels:
        for reference in kernel.root_node.walk(Reference):
            if reference.symbol.name == 'may_day_flag':
                kernels.kernels.remove(kernel)

    # disable streams
    for kernel in kernels.kernels:
        kernel.acc_async_stream = None

    # gen acc
    kernels.make_acc_tranformation(False)
    if merge_joinable:
        kernels.merge_joinable_kernels()

    # apply collaspse
    if collapse:
        apply_acc_loop_and_collapse(kernels, options)

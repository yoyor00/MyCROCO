##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

'''A transformation script that seeks to apply OpenACC KERNELS directives to
NEMO style code. In order to use it you must first install PSyclone. See
README.md in the top-level directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s <this_script> <target_source_file>

The transformation script attempts to insert Kernels directives at the
highest possible location(s) in the schedule tree (i.e. to enclose as
much code as possible in each Kernels region).

'''

from poseidon.dsl.helper import *
from poseidon.dsl.kernel import PsycloneACCKernelsDirective
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.if_block import IfBlock
from psyclone.psyir.nodes.call import Call
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.nodes import ACCDirective, RegionDirective
from psyclone.f2pygen import DirectiveGen
from poseidon_extensions import ACCWaitDirective, ACCSetDeviceNumDirective
from psyclone.transformations import ACCEnterDataTrans
from psyclone.psyGen import Transformation
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.nodes import ACCDataDirective, ACCDirective, \
    ACCEnterDataDirective, ACCKernelsDirective, ACCLoopDirective, \
    ACCParallelDirective, ACCRoutineDirective, Assignment, CodeBlock, \
    Directive, Loop, Node, OMPDeclareTargetDirective, \
    OMPDirective, OMPMasterDirective, \
    OMPParallelDirective, OMPParallelDoDirective, OMPSerialDirective, \
    OMPSingleDirective, OMPTaskloopDirective, PSyDataNode, Reference, \
    Return, Routine, Schedule, ACCUpdateDirective, IntrinsicCall
from utils import add_kernels, normalise_loops, \
    insert_explicit_loop_parallelism
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from psyclone.psyir.transformations import ACCUpdateTrans
from poseidon_extensions import CrocoACCRaiseKernelThroughIf

from psyclone.core import Signature
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir

def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    for invoke in psy.invokes.invoke_list:
        # Convert array and range syntax to explicit loops
        normalise_loops(
            invoke.schedule,
            unwrap_array_ranges=True,
            hoist_expressions=False,
        )

    kernels = extract_kernels_from_psyir(psy.container)
    kernels.make_acc_tranformation(False)
    kernels.merge_joinable_kernels()

    for kernel in kernels.kernels:
        loop = kernel.root_node
        try:
            ACCLoopTrans().apply(loop)
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

    routines = psy.container.walk(Routine)
    for routine in routines:
        if routine.name == "step2D_FB_tile":
            ifs = psy.container.walk(IfBlock)
            writer = FortranWriter ()
            for last_if in reversed(ifs):
                if writer(last_if.condition) == "iif == nfast":
                    last_if.if_body.children.insert(0, ACCWaitDirective(wait_queue=1))
                    break
        elif routine.name == "step2d":
            calls = psy.container.walk(Call)
            # to skip the IntrinsicCalls
            for call in calls:
                if not isinstance(call, IntrinsicCall):
                    break
            pos = call.position
            call.parent.children.insert(pos, ACCSetDeviceNumDirective(device_num='tile'))

        #################################################################
        vars_to_fetch = []
        for ref in invoke.schedule.walk(Reference):
            if ref.name.lower() in ['dc1d','cd1d','ffc1d','fc1d','cf1d','dr1d','dz1d','bc1d']:
                vars_to_fetch.append(ref.name.lower())
        vars_to_fetch = list(set(vars_to_fetch))

        if len(vars_to_fetch) > 0:
            # build sigs
            signatures = []
            for var in vars_to_fetch:
                print(var)
                signatures.append(Signature(var))
            print(signatures)

            ACCEnterDataTrans().apply(invoke.schedule, options={'signatures': signatures})
            for dir in invoke.schedule.walk(AccEnterDataDir):
                dir.sig_set = signatures

        #################################################################
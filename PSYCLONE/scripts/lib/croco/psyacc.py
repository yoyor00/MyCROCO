##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
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


##########################################################
# python
import os
# poesidon
from scripts.lib.poseidon.dsl.helper import extract_kernels_from_psyir
from scripts.lib.extensions import loops
# intenral
from scripts.lib.extensions.directives import ACCWaitDirective, ACCSetDeviceNumDirective
from scripts.lib.from_nemo_psyclone.utils import normalise_loops
# psyclone
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.if_block import IfBlock
from psyclone.psyir.nodes.call import Call
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.nodes import Loop, Routine, IntrinsicCall
from psyclone.transformations import ACCLoopTrans
from psyclone.psyGen import Invoke

##########################################################
def apply_psyclone_original_trans(routine: Routine, dump_snippets = False):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    # Iterate over all kernel loops and normalize them
    for invoke in routine.walk(Invoke):
        # Convert array and range syntax to explicit loops
        normalise_loops(
            invoke.schedule,
            unwrap_array_ranges=True,   # Replace [:] by concrete ranges
            hoist_expressions=False,    # TODO: comment
        )

    # Extract all kernel loops
    kernels = extract_kernels_from_psyir(routine)

    # Attach ACC parallelization information to be filled in later on
    kernels.make_acc_tranformation(False)

    # Kernel loop fusion, if same ranges, etc. and no dependencies
    kernels.merge_joinable_kernels()

    # Iterate over all kernel loops
    for kernel in kernels.kernels:
        loop = kernel.root_node

        ############################################################
        # to generate source to build unit tests
        if dump_snippets:
            loops.dump_node_as_source(loop, "kernel-openacc-origin")

        # Apply 'acc' to known kernel loops types. Skip kernel loops where this transformation can't apply (dependencies, etc.)
        try:
            ACCLoopTrans().apply(loop)
        except TransformationError as err:
            # This loop can not be transformed, proceed to next loop
            print("Loop not parallelised because:", str(err))
            continue

        # Count the number of perfectly nested loops (no if branches, etc.)
        # If perfectly nested loops, collapse them with 'acc' statement
        num_nested_loops = 0
        next_loop = loop
        while isinstance(next_loop, Loop):
            num_nested_loops += 1
            if len(next_loop.loop_body.children) > 1:
                break
            next_loop = next_loop.loop_body.children[0]

        if num_nested_loops > 1:
            # Attach information about how many loops can be collapsed
            loop.parent.parent.collapse = num_nested_loops

    # get routine walker from 'psy' ir
    if routine.name == "step2D_FB_tile":
        """
        Insert some particular ACC directives:

        Because of using "acc async" on the kernel loops, we need to add a sync step here at the end.
        """
        ifs = routine.walk(IfBlock)
        writer = FortranWriter ()
        for last_if in reversed(ifs):
            if writer(last_if.condition) == "iif == nfast":
                last_if.if_body.children.insert(0, ACCWaitDirective(wait_queue=1))
                break


    elif routine.name == "step2d":
        """
        Insert some particular ACC directives:

        No clue why this is required, it just mimics the original code
        """
        calls = routine.walk(Call)
        # to skip the IntrinsicCalls
        for call in calls:
            if not isinstance(call, IntrinsicCall):
                break
        pos = call.position
        directive = ACCSetDeviceNumDirective(device_num='tile')
        call.parent.children.insert(pos, directive)

    #################################################################

    """
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
    """

    #################################################################

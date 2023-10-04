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

VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

from extensions import acc
from extensions import scratch
from extensions import loops
from poseidon.dsl.helper import *
from psyclone.psyir.nodes.routine import Routine
from psyclone.transformations import ACCEnterDataTrans
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.nodes import Loop, Node, Reference, \
    Routine
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from psyclone.core import Signature
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir, InlinedKern

def apply_acc_loop_collapse(kernels: KernelList, options: dict) -> None:
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

def apply_acc_kernel(container: Node, collapse: bool, ignore_loops: list=[], merge_joinable: bool = False, options: dict = None) -> None:
    
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
        apply_acc_loop_collapse(kernels, options)

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

def trans(psy):
    """
    A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to Croco code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    """

    # steps
    acc.add_missing_device_vars(psy.container)
    acc.set_device_tile(psy.container)

    routines = psy.container.walk(Routine)
    for routine in routines:  
        # split the top loops on k
        if routine.name == "step3d_t_tile" or routine.name == "step3d_uv2_tile" or routine.name == "rhs3d_tile" or routine.name == "pre_step3d_tile" or routine.name == 'prsgrd_tile' or routine.name == "omega_tile":
            ############################################################
            # add scratch 3d vars
            if routine.name == "step3d_t_tile":
                scratch_3d_id = 2
                for var in ['fx', 'fe', 'work']:
                    scratch.add_3d_scratch_var(psy.container, routine, var, scratch_3d_id)
                    scratch_3d_id += 1

            # add scratch 3d vars
            if routine.name == "pre_step3d_tile":
                scratch_3d_id = 4
                for var in ['fx', 'fe', 'work']:
                    scratch.add_3d_scratch_var(psy.container, routine, var, scratch_3d_id)
                    scratch_3d_id += 1

            ############################################################
            # add scratch 1d vars
            for var in VARS_1D:
                scratch.add_1d_scratch_var(routine, var)

            ############################################################
            # First extract what we need to trans (to avoid issues as we change the tree while we run over it)
            top_loop: Loop
            loops_to_trans = []
            for top_loop in routine.walk(Loop, stop_type=Loop):
                # extract nested loop indice order
                vars=loops.extract_loop_indices_order(top_loop, exclude=['itrc'])
                print(vars[0:5])
                if vars[0:3] == ['k','j','i'] or vars[0:3] == ['j','k','i'] or vars[0:3] == ['j','i','k']:
                    loops_to_trans.append(top_loop)

            # now apply
            for top_loop in loops_to_trans:
                # extract nested loop indice order
                vars=loops.extract_loop_indices_order(top_loop, exclude=['itrc'])

                # handle 'kji' loops kind
                if vars[0:3] == ['k','j','i']:
                    if routine.name == "pre_step3d_tile":
                        #TODO might look to make work, work2 to possibly fix an issue and do not apply
                        loops.handle_kji_loop(top_loop, VARS_3D)
                    else:
                        scratch.patch_scratch_3d_arrays(top_loop, VARS_3D)
                elif vars[0:3] == ['j','k','i']:
                    loops.handle_jki_loop(top_loop, VARS_1D)
                    acc.set_private_on_loop(top_loop, 'i', ['fc1d', 'cf1d', 'dc1d', 'dZ1D', 'dR1D'])
                elif vars[0:3] == ['j','i','k']:
                    loops.handle_jik_loop(top_loop, VARS_1D, do_k_loop_fuse=True)
                    acc.set_private_on_loop(top_loop, 'i', ['fc1d', 'cf1d', 'dc1d', 'bc1d'])
            
            ############################################################
            # add scratch 3d vars
            #if routine.name == "step3d_uv1_tile":
            #    acc.set_private_on_loop(routine, 'j', ['dc'])


    ####################################################################
    joinable = False
    collapse = True
    if routine.name == "rhs3d_tile":
        joinable = False
    if routine.name == "diag_tile":
        joinable = False
    if routine.name == "pre_step3d_tile":
        collapse = True
    for routine in routines:  
        if routine.name == 'diag_tile':
            apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable, options={'independent': False})
        elif routine.name != 'set_HUV1':
            apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable)

    ####################################################################
    for routine in routines:  
        # Special handline to make like by-hand version, but can also use the uv2 standard approach
        # to be confirmed the other way is also valid as it does not genereate exact same
        # code than by hand.
        if routine.name == "step3d_uv1_tile":
            acc.set_private_on_loop(routine, 'j', ['dc'])
        if routine.name == "wvlcty_tile":
            acc.set_private_on_loop(routine, 'j', ['wrk'])
        if routine.name == "set_HUV2_tile":
            acc.set_private_on_loop(routine, 'j', ['dc','fc'])
        if routine.name == "prsgrd_tile":
            acc.set_private_on_loop(routine, 'i', ['dz1d','dr1d'])

    ##################################################################
    # automatic bench to fix the private issue
    #apply_acc_fetch_vars(psy)

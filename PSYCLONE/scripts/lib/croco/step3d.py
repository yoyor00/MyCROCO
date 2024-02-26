##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# specific CROCO vars
VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

##########################################################
# python
import os
# internal
from scripts.lib.extensions import acc
from scripts.lib.extensions import kernels
from scripts.lib.extensions import scratch
from scripts.lib.extensions import loops
# poseidon
from scripts.lib.poseidon.dsl.helper import *
from scripts.lib.extensions import loops
# psyclone
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes import Loop, Routine, Container

##########################################################
def apply_step3d_routine_trans(routine: Routine, container: Container, dump_snippets=False):
    # split the top loops on k
    if routine.name == "step3d_t_tile" or routine.name == "step3d_uv2_tile" or routine.name == "rhs3d_tile" or routine.name == "pre_step3d_tile" or routine.name == 'prsgrd_tile' or routine.name == "omega_tile":
        ############################################################
        # add scratch 3d vars
        if routine.name == "step3d_t_tile" and container:
            scratch_3d_id = 2
            for var in ['fx', 'fe', 'work']:
                scratch.add_3d_scratch_var(container, routine, var, scratch_3d_id)
                scratch_3d_id += 1

        # add scratch 3d vars
        if routine.name == "pre_step3d_tile" and container:
            scratch_3d_id = 4
            for var in ['fx', 'fe', 'work']:
                scratch.add_3d_scratch_var(container, routine, var, scratch_3d_id)
                scratch_3d_id += 1

        ############################################################
        # add scratch 1d vars
        for var in VARS_1D:
            scratch.add_1d_scratch_var(routine, var)

        ############################################################
        # to generate source to build unit tests
        if dump_snippets:
            for top_loop in routine.walk(Loop, stop_type=Loop):
                loops.dump_node_as_source(top_loop, "kernel-step3d-loop-origin")

        ############################################################
        # First extract what we need to trans (to avoid issues as we change the tree while we run over it)
        top_loop: Loop
        loops_to_trans = []
        for top_loop in routine.walk(Loop, stop_type=Loop):
            # extract nested loop indice order
            vars=loops.extract_loop_indices_order(top_loop, exclude=['itrc'])
            #print(vars[0:5])
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
                    loops.handle_kji_loop(top_loop, VARS_3D, dump_snippets=dump_snippets)
                else:
                    scratch.patch_scratch_3d_arrays(top_loop, VARS_3D)
            elif vars[0:3] == ['j','k','i']:
                loops.handle_jki_loop(top_loop, VARS_1D, dump_snippets=dump_snippets)
                acc.set_private_on_loop(top_loop, 'i', ['fc1d', 'cf1d', 'dc1d', 'dZ1D', 'dR1D'])
            elif vars[0:3] == ['j','i','k']:
                loops.handle_jik_loop(top_loop, VARS_1D, do_k_loop_fuse=True, dump_snippets=dump_snippets)
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
    if routine.name == 'diag_tile':
        kernels.apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable, options={'independent': False})
    elif routine.name != 'set_HUV1':
        kernels.apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable)

    ####################################################################
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

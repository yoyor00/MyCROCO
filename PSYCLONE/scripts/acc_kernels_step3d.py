##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

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
# specific CROCO vars
VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

##########################################################
# internal
from extensions import acc
from extensions import kernels
from extensions import scratch
from extensions import loops
# poseidon
from poseidon.dsl.helper import *
# psyclone
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes import Loop, Routine

##########################################################
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
            kernels.apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable, options={'independent': False})
        elif routine.name != 'set_HUV1':
            kernels.apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable)

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

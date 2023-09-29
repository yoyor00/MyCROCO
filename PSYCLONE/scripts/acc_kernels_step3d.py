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

from helper import add_missing_device_vars
from poseidon.dsl.helper import *
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.call import Call
from poseidon_extensions import ACCSetDeviceNumDirective
from psyclone.transformations import ACCEnterDataTrans
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.symbols import DataSymbol, REAL_TYPE, INTEGER_TYPE, BOOLEAN_TYPE, ArgumentInterface
from psyclone.psyir.nodes import Loop, Node, Reference, ArrayReference, \
    Routine, Literal, BinaryOperation, ACCLoopDirective, Schedule, IntrinsicCall
from psyclone.psyir.symbols import Symbol, ArrayType
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from psyclone.core import Signature
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir, InlinedKern
from psyclone.psyir.transformations.loop_fuse_trans import LoopFuseTrans
from psyclone.psyir.transformations.loop_swap_trans import LoopSwapTrans

def set_device_tile(psy) -> None:
    """
    Mimic existing ACC code.

    Future: This can be done automatically.
    """
    routines = psy.container.walk(Routine)
    for routine in routines:
        # set device directive
        if routine.name == "step3d_t" or routine.name == "step3d_uv2":
            """
            Some hack to insert an 'acc device_num=tile'
            """
            # add acc set device
            calls = psy.container.walk(Call)
            for call in calls:
                if not isinstance(call, IntrinsicCall):
                    break
            pos = call.position
            call.parent.children.insert(pos, ACCSetDeviceNumDirective(device_num='tile'))

def add_1d_scratch_var(routine: Routine, var_name:str) -> None:
    """
    Scratch arrays are just temporary arrays to speedup things and save memory.
    This is to mimic existing ACC code.

    Future: This can be done automatically.
    """
    # dup 3D
    sym_N = routine.symbol_table.lookup("N")
    sym_1d_name = f"{var_name}1d"
    sym_1d = DataSymbol(sym_1d_name, ArrayType(REAL_TYPE, [ArrayType.ArrayBounds(Literal("0", INTEGER_TYPE), Reference(sym_N))]))

    # reg
    routine.symbol_table.add(sym_1d)

def add_3d_scratch_var(psy, routine: Routine, var_2d_name: str, scratch_id: int) -> None:
    """
    Scratch arrays are just temporary arrays to speedup things and save memory.
    This is to mimic existing ACC code.

    Future: This can be done automatically.
    """
    # dup 3D
    sym_N = routine.symbol_table.lookup("N")
    sym_Istr = routine.symbol_table.lookup("istr")
    sym_Iend = routine.symbol_table.lookup("iend")
    sym_Jstr = routine.symbol_table.lookup("jstr")
    sym_Jend = routine.symbol_table.lookup("jend")
    sym_3d_name = f"{var_2d_name}_3d"
    # TODO : avoid this by copyging from the original var
    # currently has issues because it also copy the name of the var
    sym_3d = DataSymbol(sym_3d_name, ArrayType(REAL_TYPE, [
                                                ArrayType.ArrayBounds(
                                                    BinaryOperation.create(BinaryOperation.Operator.SUB, Reference(sym_Istr), Literal("2", INTEGER_TYPE)),
                                                    BinaryOperation.create(BinaryOperation.Operator.ADD, Reference(sym_Iend), Literal("2", INTEGER_TYPE))),
                                                ArrayType.ArrayBounds(
                                                    BinaryOperation.create(BinaryOperation.Operator.SUB, Reference(sym_Jstr), Literal("2", INTEGER_TYPE)),
                                                    BinaryOperation.create(BinaryOperation.Operator.ADD, Reference(sym_Jend), Literal("2", INTEGER_TYPE))),
                                                Reference(sym_N)
                        ]))
    #sym_3d.copy_properties(sym)
    #for shape_entry in sym.shape:
    #    sym_3d.shape.append(shape_entry.copy())
    #sym_3d.shape.append(sum_N)

    # reg
    original_arg_list = routine.symbol_table.argument_list[:]
    sym_3d.interface = ArgumentInterface(ArgumentInterface.Access.READWRITE)
    routine.symbol_table.specify_argument_list(original_arg_list + [sym_3d])
    routine.symbol_table.add(sym_3d)

    # add to call
    for call in psy.container.walk(Call):
        if call.routine.name == routine.name:
            scratch_var = call.ancestor(Routine).symbol_table.lookup("A3d")
            call.children.append(ArrayReference.create(scratch_var,[
                                                                Literal("1", INTEGER_TYPE),
                                                                Literal(str(scratch_id), INTEGER_TYPE),
                                                                Reference(Symbol("trd"))
                                                                ]))

def extract_loop_indices_order(top_loop: Loop, exclude=[]) -> list:
    """
    Extract loop variable names and return them as array.
    E.g., ['l','j','k']
    """
    vars=[]
    for inloop in top_loop.walk(Loop):
        vars.append(inloop.variable.name)   # Get loop variable
    for indice in exclude:
        if indice in vars:
            vars.remove(indice)
    return vars

def get_first_loop_on(top_loop: Loop, var: str) -> Loop:
    """
    Return psy representation of loop using 'var'
    """
    while top_loop.variable.name != var:
        top_loop = top_loop.walk(Loop)[1]
    return top_loop

def detach_and_get_childs(loop: Loop) -> list:
    """
    TODO: you know better...
    """
    # extract content
    to_detach = []
    for op in loop.loop_body:
        to_detach.append(op)

    # detach all child ops to keep the nested loops as template to copy
    ops = []
    for op in to_detach:
        ops.append(op)
        op.detach()

    # ok
    return ops


def patch_scratch_1d_arrays(top_loop: Loop) -> None:
    """
    Convert 1D scratch arrays to 2D ones to avoid race conditions if parallelizing over the outer loop
    """
    for ref in top_loop.walk(ArrayReference):
        if ref.name.lower() in VARS_1D:
            new_name = ref.name.lower()+'1d'
            ref.symbol = Symbol(new_name)
            ref.indices.pop(0)


def patch_scratch_3d_arrays(top_loop: Loop) -> None:
    """
    Convert 2D scratch arrays to 3D ones to avoid race conditions if parallelizing over the outer loop

    do k = 1, n, 1
      do j = jstr, jend, 1
        do i = MAX(istr - 1, 1), MIN(iend + 2, lm + 1), 1
          fx(i,j) = t(i,j,k,nadv,itrc) - t(i - 1,j,k,nadv,itrc)

    do k = 1, n, 1
      do j = jstr, jend, 1
        do i = MAX(istr - 1, 1), MIN(iend + 2, lm + 1), 1
          fx_3d(i,j,k) = t(i,j,k,nadv,itrc) - t(i - 1,j,k,nadv,itrc)
    """
    for ref in top_loop.walk(ArrayReference):
        if ref.name.lower() in  VARS_3D:
            new_name = ref.name.lower()+'_3d'
            ref.symbol = Symbol(new_name)
            ref.indices.append(Reference(Symbol('k')))

def handle_kji_loop(top_loop: Loop) -> None:
    """
    Swapping indices of 'kji' loops to TODO
    """
    # get top k loop
    k_loop = get_first_loop_on(top_loop, 'k')

    # check ok
    assert k_loop.variable.name == 'k'

    # get infos
    k_loop_position = k_loop.position
    parent_node = k_loop.parent

    # detach all child ops to keep the k loop as template
    ops = detach_and_get_childs(k_loop)

    # remove template k_loop
    k_loop.detach()
   
    # rebuild
    op: Node
    for i, op in enumerate(ops):
        # re inject in place of k_loop
        print(parent_node.view())
        parent_node.addchild(op, k_loop_position + i)

        # TODO before downling in ifs we need to check we didn't crossed
        # usage of k in the in-between ops
        for loop in op.walk(Loop, stop_type=Loop):
            # remember loop pos
            loop_position = loop.position
            loop_parent = loop.parent

            # detach it to put in in k_loop
            loop.detach()
            new_k_loop = k_loop.copy()
            new_k_loop.loop_body.addchild(loop)

            # put in place
            if loop_parent != None:
                loop_parent.addchild(new_k_loop, loop_position)

     # patch vars
    patch_scratch_3d_arrays(top_loop)

def handle_kji_loop_old(top_loop: Loop) -> None:
    """
    TODO: obsolete
    """
    # get top k loop
    k_loop = get_first_loop_on(top_loop, 'k')

    # check ok
    assert k_loop.variable.name == 'k'

    # get infos
    k_loop_position = k_loop.position
    parent_node = k_loop.parent

    # detach all child ops to keep the k loop as template
    ops = detach_and_get_childs(k_loop)

    # rebuild N k-loops
    new_loops = [k_loop]
    pos = k_loop_position
    for i in range(len(ops) - 1):
        new_empty_loop = k_loop.copy()
        new_loops.append(new_empty_loop)

    # fill
    for i, op in enumerate(ops):
        inner_new_loop = new_loops[i]
        inner_new_loop.loop_body.addchild(op)

    # attach loops
    for i in range(len(new_loops) - 1):
        parent_node.addchild(new_loops[i+1], pos + i + 1)

    # patch vars
    patch_scratch_3d_arrays(top_loop)


def is_loop_using_var(loop: Loop, vars: list):
    """
    Check if loop is using var as index
    """
    for ref in loop.walk(Reference):
        if ref.name in vars:
            return True
    return False


def set_private_on_loop(top_loop: Node, loop_var: str, vars:list):
    """
    Add acc private on all these loops

    TODO: if it's ACC related, use acc_ prefix for this function call.
    """
    loop: Loop
    for loop in top_loop.walk(Loop):
        if loop.variable.name == loop_var and is_loop_using_var(loop, vars):
            parent = loop.parent
            pos = loop.position
            loop_directive = ACCLoopDirective(private=vars)
            loop.detach()
            loop_directive.children[0].children.append(loop)
            parent.children.insert(pos, loop_directive)

def handle_jki_loop(top_loop: Loop) -> None:
    """
    TODO: Describe it
    """
    # remove inner j loops
    for inner_i_loop in top_loop.walk(Loop):
        if inner_i_loop.variable.name == 'i':
            parent = inner_i_loop.parent
            inner_loop_position = inner_i_loop.position
            childs = []
            for op in inner_i_loop.loop_body:
                childs.append(op)
            for op in childs:
                op.detach()
            for op in childs:
                parent.addchild(op, inner_loop_position)
            inner_i_loop.detach()
    # put i loop on top to make j,i, k 
    childs = []
    for op in top_loop.loop_body:
        childs.append(op)
    for op in childs:
        op.detach()
    for op in childs:
        inner_i_loop.loop_body.addchild(op)
    top_loop.loop_body.addchild(inner_i_loop)

    # patch arrays
    patch_scratch_1d_arrays(top_loop)

    # add private to i loops
    set_private_on_loop(top_loop, 'i', ['fc1d', 'cf1d', 'dc1d', 'dZ1D', 'dR1D'])

def handle_jik_loop(top_loop: Loop, do_k_loop_fuse: bool = True) -> None:
    """
    Describe what it does
    """
    # remove inner j loops
    #for inner_i_loop in top_loop.walk(Loop):
    #    if inner_i_loop.variable.name == 'i':
    #        parent = inner_i_loop.parent
    #        inner_loop_position = inner_i_loop.position
    #        childs = []
    #        for op in inner_i_loop.loop_body:
    #            childs.append(op)
    #        for op in childs:
    #            op.detach()
    #        for op in childs:
    #            parent.addchild(op, inner_loop_position)
    #        inner_i_loop.detach()
    ## put i loop on top to make j,i, k 
    #childs = []
    #for op in top_loop.loop_body:
    #    childs.append(op)
    #for op in childs:
    #    op.detach()
    #for op in childs:
    #    inner_i_loop.loop_body.addchild(op)
    #top_loop.loop_body.addchild(inner_i_loop)

    # swap the i inner loop
    for inner_i_loop in top_loop.walk(Loop):
        parent_loop = inner_i_loop.ancestor(Loop)
        if parent_loop != None and parent_loop.variable.name == 'k':
            if inner_i_loop.variable.name == 'i':
                LoopSwapTrans().apply(parent_loop)

    # remove inline kernel which generate bug
    for ikernel in top_loop.walk(InlinedKern):
        ikernel.parent.replace_with(Schedule(children=ikernel.children[0].pop_all_children()))

    # fuse k loops
    if do_k_loop_fuse:
        lst = []
        for inner_i_loop in top_loop.walk(Loop):
            if inner_i_loop.variable.name == 'i':
                lst.append(inner_i_loop)
        top_id = 0
        next_id = top_id + 1
        while next_id < len(lst):
            if lst[top_id].parent == lst[next_id].parent and lst[top_id].position == lst[next_id].position - 1:
                merge = True
                LoopFuseTrans().apply(lst[top_id], lst[next_id])
                next_id += 1
            else:
                top_id = next_id
                next_id += 1

    # patch arrays
    patch_scratch_1d_arrays(top_loop)

    # add private to i loops
    set_private_on_loop(top_loop, 'i', ['fc1d', 'cf1d', 'dc1d', 'bc1d'])


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
                print(var)
                signatures.append(Signature(var))
            print(signatures)

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
    add_missing_device_vars(psy.container)
    set_device_tile(psy)

    print(psy.container.view())

    routines = psy.container.walk(Routine)
    for routine in routines:  
        # split the top loops on k
        if routine.name == "step3d_t_tile" or routine.name == "step3d_uv2_tile" or routine.name == "rhs3d_tile" or routine.name == "pre_step3d_tile" or routine.name == 'prsgrd_tile' or routine.name == "omega_tile":
            ############################################################
            # add scratch 3d vars
            if routine.name == "step3d_t_tile":
                scratch_3d_id = 2
                for var in ['fx', 'fe', 'work']:
                    add_3d_scratch_var(psy, routine, var, scratch_3d_id)
                    scratch_3d_id += 1

            # add scratch 3d vars
            if routine.name == "pre_step3d_tile":
                scratch_3d_id = 4
                for var in ['fx', 'fe', 'work']:
                    add_3d_scratch_var(psy, routine, var, scratch_3d_id)
                    scratch_3d_id += 1

            ############################################################
            # add scratch 1d vars
            for var in VARS_1D:
                add_1d_scratch_var(routine, var)

            ############################################################
            # handle loop kinds
            top_loop: Loop
            for top_loop in routine.walk(Loop, stop_type=Loop):
                # extract nested loop indice order
                vars=extract_loop_indices_order(top_loop, exclude=['itrc'])
                print(vars[0:5])

                # handle 'kji' loops kind
                if vars[0:3] == ['k','j','i']:
                    if routine.name == "pre_step3d_tile":
                        #TODO might look to make work, work2 to possibly fix an issue and do not apply
                        handle_kji_loop(top_loop)
                    else:
                        patch_scratch_3d_arrays(top_loop)
                elif vars[0:3] == ['j','k','i']:
                    handle_jki_loop(top_loop)
                elif vars[0:3] == ['j','i','k']:
                    handle_jik_loop(top_loop, do_k_loop_fuse=True)
            
            ############################################################
            # add scratch 3d vars
            #if routine.name == "step3d_uv1_tile":
            #    set_private_on_loop(routine, 'j', ['dc'])


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
            set_private_on_loop(routine, 'j', ['dc'])
        if routine.name == "wvlcty_tile":
            set_private_on_loop(routine, 'j', ['wrk'])
        if routine.name == "set_HUV2_tile":
            set_private_on_loop(routine, 'j', ['dc','fc'])
        if routine.name == "prsgrd_tile":
            set_private_on_loop(routine, 'i', ['dz1d','dr1d'])

    ##################################################################
    # automatic bench to fix the private issue
    #apply_acc_fetch_vars(psy)

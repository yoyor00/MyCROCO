##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Some helper function to make change the loops in CROCO.

Remark
------
There is a lot of things currently hard coded and which requires to be made
generic intead.
'''

##########################################################
import os
from .scratch import patch_scratch_3d_arrays, patch_scratch_1d_arrays
from .acc import set_private_on_loop
from .loops_helpers import *
from psyclone.psyir.nodes import Node, Reference, Loop, Routine, Schedule
from psyclone.psyir.transformations.loop_fuse_trans import LoopFuseTrans
from psyclone.psyir.transformations.loop_swap_trans import LoopSwapTrans
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir, InlinedKern
from psyclone.psyir.backend.fortran import FortranWriter

##########################################################
ENABLE_SNIPPET_DUMPS=False

##########################################################
def dump_node_as_source(node: Node, variant: str = "origin"):
    '''
    Generated the snippets to reuse in unit test. Needs to be copied
    by hand into scripts/extension/tests/croco-loops-snippets/ in right
    directory.
    '''
    # get parent routine
    routine = node.ancestor(Routine)

    # calc ID
    childs = routine.walk(Node)
    cursor = childs[0]
    id = 0
    while cursor is not node:
        id += 1
        cursor = childs[id]

    # create directory if not exist
    path = os.path.expanduser(f'./loops-snippets/{variant}')
    os.makedirs(path, exist_ok=True)

    # open file
    fname = os.path.join(path, f"snippet-{routine.name}-{id}.F")
    with open(fname, "w+") as fp:
        source = FortranWriter()(node)
        fp.write(source)

##########################################################
def handle_kji_loop(top_loop: Loop, patch_scratch_vars: list) -> None:
    """
    Swapping indices of 'kji' loops to TODO

    - pre_step3d_tile
    """

    # to generate source to build unit tests
    if ENABLE_SNIPPET_DUMPS:
        dump_node_as_source(top_loop, "kji-origin")

    # help finding where it applies (debug)
    #if not top_loop.ancestor(Routine).name in ['pre_step3d_tile']:
    #    raise Exception(f"Applying on {top_loop.ancestor(Routine).name}")

    # patch vars
    patch_scratch_3d_arrays(top_loop, patch_scratch_vars)

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
        parent_node.addchild(op, k_loop_position + i)

        # check error, cannot go through (might need to handle it if ithappens
        # somewhere in croco). Typically if there is a "if (k....)" (or other
        # computations depending on k).
        # in this cas, we might want to down the top k loop to the if and not
        # down more into the down loop.
        for reference in op.walk(Reference, stop_type=Loop):
            if reference.symbol.name == 'k':
                raise Exception(f"Cannot down the k-loop down to the next loop "
                                f"level because it"
                                f"crosses some k-use in the path !")

        # swap down the k loop to make smaller sub independant kernels
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

##########################################################
def handle_jki_loop(top_loop: Loop, scratch_1d_vars:list) -> None:
    """
    TODO: Describe it

    - step3d_t_tile
    - acc_kernels_step3d
    - pre_step3d_tile
    - prsgrd_tile
    - rhs3d_tile
    """

    # to generate source to build unit tests
    if ENABLE_SNIPPET_DUMPS:
        dump_node_as_source(top_loop, "jki-origin")

    # help finding where it applies (debug)
    #if not top_loop.ancestor(Routine).name in ['step3d_t_tile', 'acc_kernels_step3d', 'pre_step3d_tile', 'prsgrd_tile', 'rhs3d_tile']:
    #    raise Exception(f"Applying on {top_loop.ancestor(Routine).name}")

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

##########################################################
def handle_jik_loop(top_loop: Loop, scratch_1d_vars:list, do_k_loop_fuse: bool = True) -> None:
    """
    Describe what it does

    - step3d_uv2_tile
    - rhs3d_tile
    - omega_tile
    """

    # to generate source to build unit tests
    if ENABLE_SNIPPET_DUMPS:
        dump_node_as_source(top_loop, "jik-origin")

    # patch arrays
    patch_scratch_1d_arrays(top_loop, scratch_1d_vars)

    # help finding where it applies (debug)
    #if not top_loop.ancestor(Routine).name in ['step3d_uv2_tile', 'rhs3d_tile', 'omega_tile']:
    #    raise Exception(f"Applying on {top_loop.ancestor(Routine).name}")

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

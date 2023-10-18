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
# psyclone
from psyclone.psyir.nodes import Reference, Loop

##########################################################
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

##########################################################
def get_first_loop_on(top_loop: Loop, var: str) -> Loop:
    """
    Return psy representation of loop using 'var'
    """
    while top_loop.variable.name != var:
        top_loop = top_loop.walk(Loop)[1]
    return top_loop

##########################################################
def detach_and_get_childs(loop: Loop) -> list:
    """
    Detach all the childs from the given loop and return them as a list
    (in order to let the caller injecting them at another place).
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

##########################################################
def is_loop_using_var(loop: Loop, vars: list):
    """
    Check if loop is using var as index
    """
    for ref in loop.walk(Reference):
        if ref.name in vars:
            return True
    return False

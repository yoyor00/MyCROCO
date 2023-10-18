#!/usr/bin/env python3

from psyclone.psyir.nodes.reference import Reference
from psyclone.psyir.nodes.assignment import Assignment
from psyclone.psyir.nodes.array_reference import ArrayReference
from psyclone.psyir.nodes.node import Node
from psyclone.psyir.nodes.loop import Loop
from psyclone.psyir.nodes.routine import Routine

def extract_var_ref_list(node: Node, id: int) -> list:
    result = []
    if isinstance(node, ArrayReference):
        if isinstance(node.parent, Assignment) and id == 0:
            dep = {'var': node.name, 'mode': 'W'}
        else:
            dep = {'var': node.name, 'mode': 'R'}
        result.append(dep)
    elif isinstance(node, Reference):
        if isinstance(node.parent, Assignment) and id == 0:
            dep = {'var': node.name, 'mode': 'W'}
        else:
            dep = {'var': node.name, 'mode': 'R'}
        result.append(dep)
    for id, child in enumerate(node.children):
        result += extract_var_ref_list(child, id)
    return result

def extract_var_ref_list_reduced(node: Node, id: int = 0) -> dict:
    lst = extract_var_ref_list(node, id)
    map = {}
    for var in lst:
        var_name = var['var']
        if var_name in map:
            if var['mode'] != map[var_name]:
                map[var_name] = 'RW'
        else:
            map[var_name] = var['mode']
    return map

def is_using_var(node: Node, name: str) -> bool:
    vars = extract_var_ref_list(node, 0)
    for var in vars:
        if var['var'] == name:
            return True
    return False

def get_first_loop(node: Node) -> Loop:
    if isinstance(node, Loop):
        return node
    else:
        for child in node.children:
            res = get_first_loop(child)
            if res != None:
                return res
    return None

def get_previous(node: Node, type):
    #get previous kernel to get state
    i = node.position - 1
    while i >= 0:
        if isinstance(node.parent.children[i], type):
            return node.parent.children[i]
        i -= 1
    return None

def get_previous_2(node: Node, type):
    #get previous kernel to get state
    precedors = node.parent.walk(type, stop_type=type)
    pos = -1
    for id, precedor in enumerate(precedors):
        if precedor.position < node.position and precedor.ancestor(Routine) == precedor.ancestor(Routine):
            pos = id
    if pos == -1:
        return None
    else:
        return precedors[pos]

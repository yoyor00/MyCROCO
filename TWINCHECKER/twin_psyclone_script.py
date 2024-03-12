# psyclone
from psyclone.psyir.nodes import Routine, Loop, Reference, Assignment, Call, BinaryOperation, UnaryOperation, Schedule, Literal, Node, IntrinsicCall
from psyclone.psyir.symbols import RoutineSymbol, NoType, ContainerSymbol, DataType, Symbol

def gen_call(node: Node) -> Call:
    t = "INTEGER"
    for literal in node.walk((Literal, Reference, UnaryOperation, IntrinsicCall)):
        if isinstance(literal, Reference):
            literal = literal.symbol
        if "REAL" in str(literal) or 'FLOAT' in str(literal):
            t = "REAL"
    fixable = ""
    if isinstance(node, Reference) and node.is_array:
        fixable = "_fixable"
    if "REAL" == t:
        twin_check_float = RoutineSymbol(f"twin_check_double{fixable}", NoType())
    elif "INTEGER" == t:
        twin_check_float = RoutineSymbol(f"twin_check_integer{fixable}", NoType())
    else:
        raise Exception(str(literal))
        twin_check_float = RoutineSymbol("twin_check_float", NoType())
    arguments = [node.copy()] # list of References
    call = Call.create(twin_check_float, arguments)
    return call

def gen_location(node: Node) -> dict:
    parent: Node = node
    while not isinstance(parent, Schedule):
        line_pos = parent.position
        parent = parent.parent
    return {"parent": parent, "pos": line_pos}

def trans(psy):
    routines = psy.container.walk(Routine)
    for routine in routines:
        module_symbol = ContainerSymbol("twin_checker", True)
        routine.symbol_table.add(module_symbol)
        for loop in routine.walk(Loop, stop_type=Loop):
            to_insert = []
            for assign in loop.walk(Assignment):
                # value
                value = assign.rhs
                store_in = assign.lhs

                # value of assign
                details = []
                for ref in value.walk((Reference, BinaryOperation, UnaryOperation, IntrinsicCall)):
                    loc_insert = gen_location(ref)
                    loc_insert['call'] = gen_call(ref)
                    details.append(loc_insert)

                # rever details
                details.reverse()
                to_insert += details

                # rign assign
                loc_insert = gen_location(store_in)
                loc_insert['pos'] += 1
                loc_insert['call'] = gen_call(store_in)
                to_insert.append(loc_insert)

            for id, entry in enumerate(to_insert):
                #print("----------------------")
                #print(entry['parent'])
                #print(entry['call'].view())
                entry['parent'].addchild(entry['call'], entry['pos'] + id)

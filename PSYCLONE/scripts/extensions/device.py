##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Some helper function to make some specific transform on the CROCO code.

Remark
------
There is a lot of things currently hard coded and which requires to be made
generic intead.
'''

##########################################################
from .directives.ACCSetDeviceNumDirective import ACCSetDeviceNumDirective
from psyclone.psyir.nodes import Node, Routine, Literal
from psyclone.psyir.symbols import DataSymbol, INTEGER_TYPE, BOOLEAN_TYPE
from psyclone.psyir.nodes import Call, IntrinsicCall

##########################################################
def add_missing_device_vars(root_node: Node) -> None:
    """
    Mimic existing hand-made acc versionby adding some GPU location variables
    as they did.

    Vars:
     - my_acc_device
     - compute_on_device

    Remark: 
    -------
    We might skip this in the end version as it is probably useless.

    Arguments:
    ----------
    psy: Node
        The root node of the parsed tree in intermediate representation so we can
        reshape it as a transformation pass.
    """
    routines: Routine
    routines = root_node.walk(Routine)
    for routine in routines:
        print(routine.name)
        # set device directive
        if routine.name == "step3d_t" or routine.name == "step3d_t_tile" or routine.name == "step3d_uv2" or routine.name == "step3d_uv2_tile":
            # add const integer 'my_acc_device' = 0
            symbol = DataSymbol("my_acc_device", INTEGER_TYPE, is_constant = True, initial_value = Literal("0", INTEGER_TYPE))
            routine.symbol_table.add(symbol)

            # add const boolean 'compute_on_device' = true
            symbol = DataSymbol("compute_on_device", BOOLEAN_TYPE, is_constant = True, initial_value = Literal("true", BOOLEAN_TYPE))
            routine.symbol_table.add(symbol)

##########################################################
def set_device_tile(root_node: Node) -> None:
    """
    Mimic existing ACC code.

    Future: This can be done automatically.
    """
    routines = root_node.walk(Routine)
    for routine in routines:
        # set device directive
        if routine.name == "step3d_t" or routine.name == "step3d_uv2":
            """
            Some hack to insert an 'acc device_num=tile'
            """
            # add acc set device
            calls = routine.walk(Call)
            for call in calls:
                if not isinstance(call, IntrinsicCall):
                    break
            pos = call.position
            call.parent.children.insert(pos, ACCSetDeviceNumDirective(device_num='tile'))

##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Some helper function to make change the scratch variables in CROCO.

Remark
------
There is a lot of things currently hard coded and which requires to be made
generic intead.
'''

##########################################################
# internal
from .directives.ACCSetDeviceNumDirective import ACCSetDeviceNumDirective
# psyclone
from psyclone.psyir.nodes import Node, Routine, Literal, Reference, ArrayReference, Loop
from psyclone.psyir.symbols import DataSymbol, INTEGER_TYPE, BOOLEAN_TYPE
from psyclone.psyir.nodes import Call, IntrinsicCall, BinaryOperation
from psyclone.psyir.symbols import Symbol, ArrayType, REAL_TYPE, ArgumentInterface

##########################################################
def add_1d_scratch_var(routine: Routine, var_name:str) -> None:
    """
    Scratch arrays are just temporary arrays to speedup things and save memory.
    This is to mimic existing ACC code.

    Arguments:
    ----------
    routine: Routine
        The root node to walk from.
    var_name: str
        The variable name to append (eg. 'dc') which will be concatenated to '1d'.

    Future:
    -------
    This can be done automatically.
    """
    # dup 3D
    sym_N = routine.symbol_table.lookup("N")
    sym_1d_name = f"{var_name}1d"
    sym_1d = DataSymbol(sym_1d_name, ArrayType(REAL_TYPE, [ArrayType.ArrayBounds(Literal("0", INTEGER_TYPE), Reference(sym_N))]))

    # reg
    routine.symbol_table.add(sym_1d)

##########################################################
def patch_scratch_1d_arrays(top_loop: Loop, var_names: list) -> None:
    """
    Convert 1D scratch arrays to 2D ones to avoid race conditions if parallelizing over the outer loop

    Note:
    -----
    Need to first call add_1d_scratch_var().
    """
    for ref in top_loop.walk(ArrayReference):
        if ref.name.lower() in var_names:
            new_name = ref.name.lower()+'1d'
            ref.symbol = Symbol(new_name)
            ref.indices.pop(0)

##########################################################
def add_3d_scratch_var(root_node: Node, routine: Routine, var_2d_name: str, scratch_id: int) -> None:
    """
    Scratch arrays are just temporary arrays to speedup things and save memory.
    This is to mimic existing ACC code.

    Future:
    -------
    This can be done automatically.
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
    for call in root_node.walk(Call):
        if call.routine.name == routine.name:
            scratch_var = call.ancestor(Routine).symbol_table.lookup("A3d")
            call.children.append(ArrayReference.create(scratch_var,[
                                                                Literal("1", INTEGER_TYPE),
                                                                Literal(str(scratch_id), INTEGER_TYPE),
                                                                Reference(Symbol("trd"))
                                                                ]))

##########################################################
def patch_scratch_3d_arrays(top_loop: Loop, var_names: list) -> None:
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
        if ref.name.lower() in  var_names:
            new_name = ref.name.lower()+'_3d'
            ref.symbol = Symbol(new_name)
            ref.indices.append(Reference(Symbol('k')))

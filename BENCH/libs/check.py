#!/bin/env python3
##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# imports
from netCDF4 import Dataset
import math
import numpy

##########################################################
def recurse_compare_current_dim(varname:str, shape_ref: tuple, shape_actual: tuple, cusor_ref, cursor_actual, dim_id: int, current_dims: list = []) -> None:
    '''
    Recursively compare the values of each dimensions of the mesh.

    It is usefull mostly if we want to see which value is incorrect opposite to numpy who currently say only
    what is the mistake. It also offer a different approach to crosscheck we are right with the way of using
    numpy on top of netcdf.

    Parameters:
    -----------
    varname: str
        Name of the variable we are looping over.
    shape_ref: tuple
        Size of each dimensions in ref file.
    shape_actual: tuple
        Size of each dimensions in ref file.
    cursor_ref:
        Dimension we are looping over. At first call it is a netCDF4.Variable, after if is an array in principle.
    cursor_actual:
        Dimension we are looping over. At first call it is a netCDF4.Variable, after if is an array in principle.
    dim_id: int
        The current dimension we are looping in.
    '''
    # check
    if len(shape_ref) != len(shape_actual):
        f"Shape is not same in both file for '{varname}'"

    # nothing to check
    if len(shape_ref) == 0:
        return

    # check
    assert dim_id < len(shape_ref)

    # check size of the dim
    if shape_actual[dim_id] != shape_ref[dim_id]:
        raise Exception(f"Shape is not same in both file for '{varname}'")

    # is last level
    if dim_id == len(shape_ref) - 1:
        # loop
        for i in range(shape_actual[dim_id]):
            if not math.isclose(cusor_ref[i], cursor_actual[i]):
                raise Exception(f"Non close equality in variable '{varname}': {cusor_ref[i]} !~= {cursor_actual[i]} at ({','.join(current_dims)})")
            if not cusor_ref[i] == cursor_actual[i]:
                raise Exception(f"Non strict equality in variable '{varname}' : {cusor_ref[i]} != {cursor_actual[i]} at ({','.join(current_dims)})")
    else:
        for i in range(shape_actual[dim_id]):
            recurse_compare_current_dim(varname, shape_ref, shape_actual, cusor_ref[i], cursor_actual[i], dim_id+1, current_dims=current_dims+[str(i)])

##########################################################
def compare_netcdf_variables(ref: Dataset, actual: Dataset) -> None:
    # loop on vars to check
    for var in ref.variables.keys():
        # print
        #print(f"---------------- Checking {var} -------------------")

        # extract shapes
        shape_ref = ref.variables[var].shape
        shape_actual = actual.variables[var].shape

        # check has same shape
        if shape_ref != shape_actual:
            raise Exception(f"Shape is not the same for variable '{var}'")

        # skip strings
        if var == "spherical":
            continue

        # faster way
        np_ref = numpy.array(ref.variables[var])
        np_actual = numpy.array(actual.variables[var])
        need_value_compare = False
        if not numpy.allclose(np_ref, np_actual):
            need_value_compare = True
        if (np_ref != np_actual).any():
            need_value_compare = True

        # ------------ if needed for debug
        # Note : kept if we want to see the exact failing value for debug
        # Same but by hand recusion
        if need_value_compare:
            recurse_compare_current_dim(var, shape_ref, shape_actual, ref.variables[var], actual.variables[var], 0)

            # in case it is not compared the same way we should not let go
            raise Exception(f"Seen error when comparing values for variable '{var}' ! This error message should never be reached !")

##########################################################
def compare_netcdf_files(ref_file: str, actual_file: str) -> None:
    # load
    ref = Dataset(ref_file, "r")
    actual = Dataset(actual_file, "r")

    # check & forward exception if has one
    try:
        compare_netcdf_variables(ref, actual)
    finally:
        ref.close()
        actual.close()

##########################################################
def compare_netcdf_files_old(ref_file: str, actuel_file: str):
    # load
    ref = Dataset(ref_file, "r")
    actual = Dataset(actuel_file, "r")

    # loop on vars to check
    for var in ref.variables.keys():
        # print
        #print(f"---------------- Checking {var} -------------------")

        # extract shapes
        shape_ref = ref.variables[var].shape
        shape_actual = actual.variables[var].shape

        # check
        #print("----------------------------------------------------")
        #print(shape_ref)
        #print(shape_actual)
        #print("----------------------------------------------------")
        assert shape_ref == shape_actual

        # compare one var
        if len(shape_ref) == 3:
            for i in range(shape_actual[0]):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                for j in range(shape_actual[1]):
                    v1_j = v1_i[j]
                    v2_j = v2_i[j]
                    for k in range(shape_actual[2]):
                        v1 = v1_j[k]
                        v2 = v2_j[k]
                        assert math.isclose(v1, v2), f"({i}, {j}, {k}) => {v1} != {v2}"
                        assert v1 == v2, f"({i}, {j}, {k}) => {v1} != {v2}"
        elif len(shape_ref) == 2:
            for i in range(shape_actual[0]):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                for j in range(shape_actual[1]):
                    v1_j = v1_i[j]
                    v2_j = v2_i[j]
                    assert math.isclose(v1_j, v2_j), f"({i}, {j}) => {v1_j} != {v2_j}"
                    assert v1_j == v2_j, f"({i}, {j}, {k}) => {v1} != {v2}"
        elif len(shape_ref) == 1:
            for i in range(shape_actual[0]):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                assert math.isclose(v1_i, v2_i), f"({i}, {j}) => {v1_i} != {v2_i}"
                assert v1_i == v2_i, f"({i}, {j}, {k}) => {v1} != {v2}"

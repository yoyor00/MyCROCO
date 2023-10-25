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

##########################################################
def compare_netcdf_files(ref_file: str, actuel_file: str):
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
            for i in range(min(shape_actual[0], shape_ref[0])):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                for j in range(min(shape_actual[1], shape_ref[1])):
                    v1_j = v1_i[j]
                    v2_j = v2_i[j]
                    for k in range(min(shape_actual[2], shape_ref[2])):
                        v1 = v2_j[k]
                        v2 = v2_j[k]
                        assert math.isclose(v1, v2), f"({i}, {j}, {k}) => {v1} != {v2}"
                        assert v1 == v2, f"({i}, {j}, {k}) => {v1} != {v2}"
        elif len(shape_ref) == 2:
            for i in range(min(shape_actual[0], shape_ref[0])):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                for j in range(min(shape_actual[1], shape_ref[1])):
                    v1_j = v1_i[j]
                    v2_j = v2_i[j]
                    assert math.isclose(v1_j, v2_j), f"({i}, {j}) => {v1_j} != {v2_j}"
                    assert v1_j == v2_j, f"({i}, {j}, {k}) => {v1} != {v2}"
        elif len(shape_ref) == 1:
            for i in range(min(shape_actual[0], shape_ref[0])):
                v1_i = ref.variables[var][i]
                v2_i = actual.variables[var][i]
                assert math.isclose(v1_i, v2_i), f"({i}, {j}) => {v1_i} != {v2_i}"
                assert v1_i == v2_i, f"({i}, {j}, {k}) => {v1} != {v2}"

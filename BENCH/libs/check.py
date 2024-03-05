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
class CompareErrorLogger:
    def __init__(self, max_stored: int = 10, max_total: int = 50):
        self.count = 0
        self.var_error_logs = {}
        self.var_error_count = {}
        self.max_stored = max_stored
        self.max_total = max_total

    def append(self, varname: str, ref, actual, coord: list, is_strict_compare: bool = True):
        # count it
        self.count += 1

        # first seen for this var
        if not varname in self.var_error_logs:
            self.var_error_logs[varname] = []
            self.var_error_count[varname] = 0

        # count & get log
        self.var_error_count[varname] += 1
        var_error_log = self.var_error_logs[varname]

        # to help readability with right operator displayed
        if is_strict_compare:
            compare_name = "strict"
            operator = "!="
        else:
            compare_name = "close"
            operator = "!~="

        # calc diff
        diff = abs(ref - actual)

        # extract some meaning on limits
        can_still_log_var = (len(var_error_log) <= self.max_stored)

        # if can still log
        if can_still_log_var:
            var_error_log.append(f"Non {compare_name} equality in variable '{varname}' at ({','.join(coord)}): ref {operator} actual : {ref} {operator} {actual} (diff={diff})")

    def has_error(self) -> bool:
        return self.count > 0

    def __str__(self):
        # prepare some vars
        error_var_names = ', '.join(self.var_error_logs.keys())
        total_error_count = self.count

        # prep per variable messages
        var_messages = []
        log_count = 0
        for varname, log in self.var_error_logs.items():
            var_messages.append(f"-------------------- {varname} ---------------------")
            for entry in log:
                log_count += 1
                if log_count < self.max_total:
                    var_messages.append(entry)
                elif log_count == self.max_total:
                    var_messages.append(f".................... too many errors (>{self.max_total}), stop logging details .............")

        # assemble details
        details = '\n'.join(var_messages)

        # build full message
        message = f"Found {total_error_count} errors in : {error_var_names}\n{details}"
        
        # ok
        return message

##########################################################
def recurse_compare_current_dim(error_log: CompareErrorLogger, varname:str, shape_ref: tuple, shape_actual: tuple, cusor_ref, cursor_actual, dim_id: int, current_dims: list = []) -> None:
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
                error_log.append(varname, cusor_ref[i], cursor_actual[i], current_dims+[str(i)])
            if not cusor_ref[i] == cursor_actual[i]:
                error_log.append(varname, cusor_ref[i], cursor_actual[i], current_dims+[str(i)])
    else:
        for i in range(shape_actual[dim_id]):
            recurse_compare_current_dim(error_log, varname, shape_ref, shape_actual, cusor_ref[i], cursor_actual[i], dim_id+1, current_dims=current_dims+[str(i)])

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
            # error logger
            error_logger = CompareErrorLogger(max_stored=10, max_total=50)

            # log all errors
            recurse_compare_current_dim(error_logger, var, shape_ref, shape_actual, ref.variables[var], actual.variables[var], 0)

            # in case it is not compared the same way we should not let go
            if not error_logger.has_error():
                raise Exception(f"Seen error when comparing values for variable '{var}' ! This error message should never be reached !")
            
            # log errors
            raise Exception(f"Detect some errors : \n{error_logger}")

##########################################################
def compare_netcdf_files(ref_file: str, actual_file: str) -> None:
    # load
    ref = Dataset(ref_file, "r")
    actual = Dataset(actual_file, "r")

    # check & forward exception if has one
    try:
        compare_netcdf_variables(ref, actual)
    except Exception as e:
        raise Exception(f"Error while checking\n - refere : {ref_file}\n - actual : {actual_file}\n-----------------------------------------------\n" + str(e))
    finally:
        ref.close()
        actual.close()

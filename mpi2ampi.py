#!/usr/bin/env python3
# Usage: mpi2ampi <some src>.f <output>.f

import sys
import shutil
import logging as log
from fparser import api
from fparser.common.readfortran import FortranFileReader
from fparser.two.Fortran2003 import Subroutine_Stmt, Call_Stmt, Name, \
    Actual_Arg_Spec_List, Return_Stmt
from fparser.two.parser import ParserFactory
from fparser.two.utils import walk
from fparser.two import utils

# check the presence of MPI_Irecv or MPI_Recv calls
def check_irecv(ast):
    found_mpi = False
    found_irecv = False
    found_recv = False
    for node in walk(ast):
        if isinstance(node, Call_Stmt):
            call_fun = node.items[0].string.lower()
            if call_fun.startswith('mpi'):
                found_mpi = True
            if call_fun == "mpi_irecv":
                found_irecv = True
            if call_fun == "mpi_recv":
                found_recv = True
    if (found_recv and found_irecv):
        print("mix of mpi_irecv and mpi_recv! stop.")
        sys.exit(1)
    else:
        return [found_mpi,found_irecv, found_recv]
        
def str_replace_start(start, r, s):
    return r + s[len(start):]  if s.lower().startswith(start.lower()) else s

# transform fun
def mpi_to_ampi_call(ast, arg_translat, tag_add):
    need_turn_calls = ["mpi_isend", "mpi_irecv"]
    need_turn = False
    need_turn_buffs = []
    current_subroutine = None
    returns = {}
    log.debug(type(ast))
    for node in walk(ast):
        log.debug("enter:{}".format(type(node)))
        if isinstance(node, Subroutine_Stmt):
            current_subroutine = str(node.get_name())
        if isinstance(node, Call_Stmt):
            call_fun = node.items[0].string.lower()
            
            if call_fun.startswith('mpi'):

                if call_fun in need_turn_calls:
                    need_turn = True
                    need_turn_buffs.append(str(node.items[1].items[0]))

                node.items[0].string = str_replace_start('MPI', 'AMPI',
                                                         node.items[0].string)
                litems = list(node.items[1].items)
                str_items = [str(i) for i in list(node.items[1].items)]
                for key in arg_translat:
                    if key in str_items:
                        index = str_items.index(key)
                        litems[index] = Name(arg_translat[key])
                if node.items[0].string in tag_add:
                    # AMPI direction tag comes before communicator spec
                    index = str_items.index("MPI_COMM_WORLD")
                    litems.insert(index, Name(tag_add[node.items[0].string]))
                node.items[1].items = tuple(litems)
        if isinstance(node, Return_Stmt):
            assert current_subroutine is not None
            assert not (current_subroutine in returns)
            returns[current_subroutine] = True
            #print (need_turn_buffs)
            need_turn = False
            need_turn_buffs = []
        
    
    return ast

if __name__ == '__main__':

    # debug used in fparser2
    #log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    log.basicConfig(format="%(levelname)s: %(message)s")
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # read + parse
    reader = FortranFileReader(input_file)
    parser = ParserFactory().create() # 2003
    ast = parser(reader)

    # irecv or recv cannot be mixed in the same src
    [found_mpi, irecv_mode, recv_mode] = check_irecv(ast)

    if True:
        # MPI -> AMPI translation
        arg_translat = { "MPI_DOUBLE_PRECISION": "AMPI_ADOUBLE_PRECISION" } #etc.
        
        ampi_tag_irecv = { "AMPI_Send": "AMPI_TO_IRECV_WAIT",
                           "AMPI_Irecv": "AMPI_FROM_SEND" }
        ampi_tag_recv = { "AMPI_Send": "AMPI_TO_RECV",
                          "AMPI_Recv": "AMPI_FROM_SEND" }

        # etc.

        if(irecv_mode):
            ast = mpi_to_ampi_call(ast, arg_translat, ampi_tag_irecv)
            
        elif(recv_mode):
            ast = mpi_to_ampi_call(ast, arg_translat, ampi_tag_recv)

        # output
        try:
            with open(output_file, "w") as f:
                f.write(ast.tofortran(tab='      ',isfix=True))
        except:
            print("cannot write to {}".format(output_file))
            sys.exit(1)

    else:
        shutil.copy(input_file, output_file)


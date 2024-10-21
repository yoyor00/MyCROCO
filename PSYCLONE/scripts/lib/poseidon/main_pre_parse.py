#!/usr/bin/env python3

import sys
import contextlib
import subprocess
from argparse import ArgumentParser
from poseidon.dsl.helper import *
from poseidon.dsl.kernel_pre_scheduler import KernelPreScheduler
from poseidon.parsing.graph_generator import GraphGenerator
from psyclone.psyir.backend.fortran import FortranWriter
from deepdiff import DeepDiff
from copy import deepcopy

def main():
    # common base options
    parser = ArgumentParser()
    parser.add_argument('-o', '--output', default='-', help="Write in given output file")
    parser.add_argument('-s', '--split', action='store_true', help="Enable loop splitting")
    parser.add_argument("source_files", nargs='+', help='List of file to parse')

    # parse
    args = parser.parse_args()

    # get files
    files = args.source_files

    # build
    sched = KernelPreScheduler()
    for file in files:
        sched.add_file(file)

    # reshape loops
    if args.split:
        sched.make_loop_splitting_assignements()

    # schedule
    sched.calculate2()

    # save
    if args.output == '-':
        sched.print()
    else:
        sched.save(args.output)

    # apply
    sched.apply()

    # try
    print("================================================")
    for file in files:
        # gen
        sched.files[file]['kernels'].apply_sched_acc_async_wait()
        #sched.files[file]['kernels'].calc_acc_async_wait(True)
        sched.files[file]['kernels'].make_acc_tranformation(False, async_deps=True)
        res = FortranWriter()(sched.files[file]['psyir'])

        print("-------------------------------------------------")

        sched2 = KernelPreScheduler()
        sched2.add_file(file)
        if args.split:
            sched2.make_loop_splitting_assignements()
        sched2.load(args.output)
        sched2.apply()

        print("-------------------------------------------------")

        sched2.files[file]['kernels'].apply_sched_acc_async_wait()
        #sched2.files[file]['kernels'].calc_acc_async_wait(True)
        sched2.files[file]['kernels'].make_acc_tranformation(False, async_deps=True)
        res2 = FortranWriter()(sched2.files[file]['psyir'])

        res_lines = res.split("\n")
        res2_lines = res2.split("\n")
        for i, val in enumerate(res_lines):
            if res_lines[i] != res2_lines[i]:
                print(f" - {res_lines[i]}")
                print(f" + {res2_lines[i]}")

        assert res2 == res


    # dump
    #sched.print_codes()
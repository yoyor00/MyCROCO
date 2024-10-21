#!/usr/bin/env python3

import sys
import contextlib
import subprocess
from argparse import ArgumentParser
from poseidon.dsl.helper import *
from poseidon.parsing.graph_generator import GraphGenerator
from psyclone.psyir.backend.fortran import FortranWriter

def parse_program_arguments() -> dict:
    # common base options
    parser = ArgumentParser()
    parser.add_argument('-g', '--global')

    # crete subcommand tree
    subparsers = parser.add_subparsers(dest="subcommand", required=True) # this line changed

    # sub command
    kernel_summary_parser = subparsers.add_parser('kernels-summary', help="Display the extracted kernels with extra informations.")
    kernel_summary_parser.add_argument('-o', '--output', default='-', help="Output file in which to write.")
    kernel_summary_parser.add_argument('-a', '--acc', action='store_true', help="If enable ACC kernel insertion.")
    kernel_summary_parser.add_argument('-s', '--split', action='store_true', help="If enable kernel loop splitting.")
    kernel_summary_parser.add_argument('SOURCE_FILE', nargs=1)

    # sub command
    kernel_deps_parser = subparsers.add_parser('kernels-deps', help="Generate the variables dependency diagram.")
    kernel_deps_parser.add_argument('-o', '--output', default='-', help="Output file in which to write.")
    kernel_deps_parser.add_argument('-f', '--format', default='dot', help="Output format to generate.")
    kernel_deps_parser.add_argument('SOURCE_FILE', nargs=1)

    # sub command
    kernel_steps_parser = subparsers.add_parser('kernels-steps', help="Generate the variables dependency diagram.")
    kernel_steps_parser.add_argument('-o', '--output', default='-', help="Output file in which to write.")
    kernel_steps_parser.add_argument('-f', '--format', default='dot', help="Output format to generate.")
    kernel_steps_parser.add_argument('SOURCE_FILE', nargs=1)

    # sub command
    all_parser = subparsers.add_parser('all', help="Generate all the kernel rendering in one go.")
    all_parser.add_argument('SOURCE_FILE', nargs=1)

    # sub command
    tree_parser = subparsers.add_parser('tree', help="Display the instruction tree.")
    tree_parser.add_argument('-o', '--output', default='-', help="Output file in which to write.")
    tree_parser.add_argument('SOURCE_FILE', nargs=1)

    # sub command
    regen_parser = subparsers.add_parser('regen', help="Regenerate the code.")
    regen_parser.add_argument('-o', '--output', default='-', help="Output file in which to write.")
    regen_parser.add_argument('-a', '--acc', action='store_true', help="If enable ACC kernel insertion.")
    regen_parser.add_argument('-s', '--split', action='store_true', help="If enable kernel loop splitting.")
    regen_parser.add_argument('-f', '--fuse', action='store_true', help="If enable kernel loop fusion.")
    regen_parser.add_argument('SOURCE_FILE', nargs=1)

    # parse
    return parser.parse_args()

# From : https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
@contextlib.contextmanager
def smart_open(filename, mode):
    if filename and filename != '-':
        fh = open(filename, mode)
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def render_dot_as_format(filename: str, dot: str, format: str):
    with smart_open(filename, 'w+') as fp:
        if format == 'dot':
            fp.write(dot)
        else:
            proc = subprocess.Popen(["dot", f"-T{format}"], stdout=fp, stdin=subprocess.PIPE)
            proc.stdin.write(dot.encode())
            proc.stdin.close()
            proc.wait()

def dump_kernels_summary(input_file: str, output_file: str = '-', acc: bool = False, split: bool = False):
    kernels = extract_kernels_from_file(input_file)

    if acc:
        kernels.make_acc_tranformation(split)
    if split:
        assert acc == True
        kernels.make_loop_splitting()

    with smart_open(output_file, 'w+') as fp:
        fp.write(kernels.render_summary())

def dump_kernels_deps(input_file: str, output_file: str = '-', format: str = 'dot'):
    kernels = extract_kernels_from_file(input_file)
    render_dot_as_format(output_file, kernels.render_dep_graph(), format)

def dump_kernels_steps(input_file: str, output_file: str = '-', format: str = 'dot'):
    kernels = extract_kernels_from_file(input_file)
    render_dot_as_format(output_file, kernels.render_step_graph(single_link=True), format)

def dump_kernels_all(input_file: str):
    kernels = extract_kernels_from_file(input_file)
    print(f" - Generating {input_file}.deps.dot ...")
    render_dot_as_format(f"{input_file}.deps.dot", kernels.render_dep_graph(), 'dot')
    print(f" - Generating {input_file}.deps.dot.svg ...")
    render_dot_as_format(f"{input_file}.deps.dot.svg", kernels.render_dep_graph(), 'svg')
    print(f" - Generating {input_file}.deps.dot.png ...")
    render_dot_as_format(f"{input_file}.deps.dot.png", kernels.render_dep_graph(), 'png')

    print(f" - Generating {input_file}.kernel.deps.dot ...")
    render_dot_as_format(f"{input_file}.kernel.deps.dot", kernels.render_dep_graph_kernels(), 'dot')
    print(f" - Generating {input_file}.kernel.deps.dot.svg ...")
    render_dot_as_format(f"{input_file}.kernel.deps.dot.svg", kernels.render_dep_graph_kernels(), 'svg')
    print(f" - Generating {input_file}.kernel.deps.dot.png ...")
    render_dot_as_format(f"{input_file}.kernel.deps.dot.png", kernels.render_dep_graph_kernels(), 'png')

    print(f" - Generating {input_file}.kernel.order.deps.dot ...")
    render_dot_as_format(f"{input_file}.kernel.order.deps.dot", kernels.render_dep_graph_kernels(exec_order_links=True), 'dot')
    print(f" - Generating {input_file}.kernel.order.deps.dot.svg ...")
    render_dot_as_format(f"{input_file}.kernel.order.deps.dot.svg", kernels.render_dep_graph_kernels(exec_order_links=True), 'svg')
    print(f" - Generating {input_file}.kernel.order.deps.dot.png ...")
    render_dot_as_format(f"{input_file}.kernel.order.deps.dot.png", kernels.render_dep_graph_kernels(exec_order_links=True), 'png')

    print(f" - Generating {input_file}.vars.deps.dot ...")
    render_dot_as_format(f"{input_file}.vars.deps.dot", kernels.render_dep_graph_vars(), 'dot')
    print(f" - Generating {input_file}.vars.deps.dot.svg ...")
    render_dot_as_format(f"{input_file}.vars.deps.dot.svg", kernels.render_dep_graph_vars(), 'svg')
    print(f" - Generating {input_file}.vars.deps.dot.png ...")
    render_dot_as_format(f"{input_file}.vars.deps.dot.png", kernels.render_dep_graph_vars(), 'png')

    print(f" - Generating {input_file}.steps.dot ...")
    render_dot_as_format(f"{input_file}.steps.dot", kernels.render_step_graph(), 'dot')
    print(f" - Generating {input_file}.steps.dot.svg ...")
    render_dot_as_format(f"{input_file}.steps.dot.svg", kernels.render_step_graph(), 'svg')
    print(f" - Generating {input_file}.steps.dot.png ...")
    render_dot_as_format(f"{input_file}.steps.dot.png", kernels.render_step_graph(), 'png')

    print(f" - Generating {input_file}.calltree.dot ...")
    render_dot_as_format(f"{input_file}.calltree.dot", gen_call_tree(input_file), 'dot')
    print(f" - Generating {input_file}.calltree.dot.svg ...")
    render_dot_as_format(f"{input_file}.calltree.dot.svg", gen_call_tree(input_file), 'svg')
    #print(f" - Generating {input_file}.calltree.dot.png ...")
    #render_dot_as_format(f"{input_file}.calltree.dot.png", gen_call_tree(input_file), 'png')

    print(f" - Generating {input_file}.summary.txt ...")
    with smart_open(f"{input_file}.summary.txt", 'w+') as fp:
        fp.write(kernels.render_summary())

def gen_call_tree(input_file: str, output_file: str = "-"):
    reader = FortranReader()
    psyir_tree = reader.psyir_from_file(input_file, free_form=False)
    gen = GraphGenerator()
    walker = Walker(gen)
    walker.walk(psyir_tree)
    return gen.graph.render()

def dump_call_tree(input_file: str, output_file: str = "-"):
    with smart_open(output_file, 'w+') as fp:
        fp.write(gen_call_tree(input_file))

def dump_regen_code(input_file: str, output_file: str = "-", acc: bool = False, split: bool = False, fuse: bool = False):
    reader = FortranReader()
    psyir_tree = reader.psyir_from_file(input_file, free_form=False)
    kernels = extract_kernels_from_psyir(psyir_tree)

    if acc:
        kernels.make_acc_tranformation(split)
        if not split:
            kernels.merge_joinable_kernels()
    if split:
        assert acc == True
        kernels.make_loop_splitting()
    if fuse:
        assert split == True
        kernels.make_loop_fusing()

    writer = FortranWriter()
    result = writer(psyir_tree)
    with smart_open(output_file, 'wb') as fp:
        fp.write(result)

def main():
    args = parse_program_arguments()

    if args.subcommand == 'kernels-summary':
        dump_kernels_summary(args.SOURCE_FILE[0], output_file=args.output, acc=args.acc, split=args.split)
    elif args.subcommand == 'kernels-deps':
        dump_kernels_deps(args.SOURCE_FILE[0], output_file=args.output, format=args.format)
    elif args.subcommand == 'kernels-steps':
        dump_kernels_steps(args.SOURCE_FILE[0], output_file=args.output, format=args.format)
    elif args.subcommand == 'all':
        dump_kernels_all(args.SOURCE_FILE[0])
    elif args.subcommand == 'tree':
        dump_call_tree(args.SOURCE_FILE[0])
    elif args.subcommand == 'regen':
        dump_regen_code(args.SOURCE_FILE[0], output_file=args.output, acc=args.acc, split=args.split, fuse=args.fuse)

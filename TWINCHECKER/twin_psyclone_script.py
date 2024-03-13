#!/usr/bin/env python3

# python
import os
import json5
import shutil
import random
import argparse
# psyclone
from psyclone.psyir.nodes import Routine, Loop, Reference, Assignment, Call, BinaryOperation, UnaryOperation, Schedule, Literal, Node, IntrinsicCall, Container, IfBlock, Operation
from psyclone.psyir.symbols import RoutineSymbol, NoType, ContainerSymbol, DataType, Symbol, ArrayType, CHARACTER_TYPE, INTEGER4_TYPE, INTEGER_UNDEF_TYPE, INTEGER8_TYPE, ScalarType
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.line_length import FortLineLength

class TwinCheckerInstrumenter:
    def __init__(self, source_file: str, dest_file: str, config: dict):
        self.fortran_writer = FortranWriter()
        self.source_file = source_file
        self.dest_file = dest_file
        self.site_ids = []
        self.config = config

        # extract some configs
        self.enable_profile_in_values = self.config.get('profile_in_values', True)
        self.enable_profile_out_values = self.config.get('profile_out_values', True)
        self.enable_profile_integers = self.config.get('profile_integers', False)

    def gen_call(self, node: Node) -> Call:
        # vars
        enable_profile_integers = self.enable_profile_integers

        t = "INTEGER"
        for literal in node.walk((Literal, Reference, UnaryOperation, IntrinsicCall)):
            if isinstance(literal, Reference):
                literal = literal.symbol
            if "REAL" in str(literal) or 'FLOAT' in str(literal):
                t = "REAL"
            if "BOOLEAN" in str(literal):
                t = "BOOLEAN"
        for literal in node.walk((Literal, Reference, UnaryOperation, IntrinsicCall)):
            if 'NINT' in str(literal):
                t = 'INTEGER'

        if t == "BOOLEAN":
            return None
        if t == "INTEGER" and not enable_profile_integers:
            return None

        fixable = ""
        if isinstance(node, Reference) and node.is_array:
            fixable = "_fixable"

        if "REAL" == t:
            twin_check_float = RoutineSymbol(f"twin_check_double{fixable}", NoType())
        elif "INTEGER" == t:
            twin_check_float = RoutineSymbol(f"twin_check_integer{fixable}", NoType())
        elif "BOOLEAN" == t:
            twin_check_float = RoutineSymbol(f"twin_check_bool{fixable}", NoType())
        else:
            raise Exception(str(literal))
            twin_check_float = RoutineSymbol("twin_check_float", NoType())
        equation_str = node.debug_string()
        if 'zob' == equation_str:
            return None
        if ':' in equation_str:
            return None
        site_id = random.randint(0, 10000000000)
        self.site_ids.append(site_id)
        arguments = [
            node.copy(), 
            Literal(equation_str, ArrayType(CHARACTER_TYPE,[len(equation_str)+1])),
            Literal(str(len(equation_str)), INTEGER8_TYPE),
            Literal(str(site_id), INTEGER8_TYPE),
            Literal("__LINE__", INTEGER_UNDEF_TYPE)
        ]
        call = Call.create(twin_check_float, arguments)
        return call

    @staticmethod
    def gen_location(node: Node) -> dict:
        parent: Node = node
        while isinstance(parent, (Reference, Operation, Call, Assignment)):
            line_pos = parent.position
            parent = parent.parent
        return {"parent": parent, "pos": line_pos}
    
    def apply_on_node(self, node: Node, stop_type = None):
        # some config needs
        enable_profile_in_values = self.enable_profile_in_values
        enable_profile_out_values = self.enable_profile_out_values

        # vars
        alread_checked = []
        to_insert = []
        loop = node

        # search assigns
        for assign in loop.walk(Assignment, stop_type=stop_type):
            # value
            value = assign.rhs
            store_in = assign.lhs

            # value of assign
            if enable_profile_in_values:
                details = []
                ref: Node
                for ref in value.walk((Reference, BinaryOperation, UnaryOperation, IntrinsicCall)):
                    as_str = ref.debug_string()
                    if as_str in alread_checked:
                        continue
                    loc_insert = self.gen_location(ref)
                    loc_insert['call'] = self.gen_call(ref)
                    details.append(loc_insert)
                    alread_checked.append(as_str)

                # rever details
                details.reverse()

                # merge
                to_insert += details

            # out assign
            if enable_profile_out_values:
                loc_insert = self.gen_location(store_in)
                loc_insert['pos'] += 1
                loc_insert['call'] = self.gen_call(store_in)
                to_insert.append(loc_insert)

        remember = {}
        for entry in to_insert:
            #print("----------------------")
            #print(entry['parent'])
            #print(entry['call'].debug_view())
            #print("----------------------")
            if entry['call'] != None:
                parent: Node = entry['parent']
                addr = parent.__repr__()
                offset = remember.get(addr, 0)
                parent.addchild(entry['call'], index=entry['pos'] + offset)
                offset += 1
                remember[addr] = offset

    def apply(self, container: Container):
        # vars
        dest_file = self.dest_file.replace('/home/svalat/Projects/minicroco/', './') #@todo: fix to keep full name but as isusue with line reshapingm
        routines = container.walk(Routine)

        for routine in routines:
            module_symbol = ContainerSymbol("twin_checker", True)
            routine.symbol_table.add(module_symbol)

            self.apply_on_node(routine, stop_type=Loop)
            for loop in routine.walk(Loop, stop_type=Loop):
                self.apply_on_node(loop)

            for id in self.site_ids:
                call = Call.create(RoutineSymbol(f"twin_register_site", NoType()), [
                    Literal(str(id), INTEGER8_TYPE),
                    Literal(dest_file, ArrayType(CHARACTER_TYPE,[len(dest_file)+1])),
                    Literal(str(len(dest_file)), INTEGER8_TYPE),
                ])
                #print(container.view())
                routine.addchild(call, 0)

def trans(psy):
    instumenter = TwinCheckerInstrumenter()
    instumenter.apply(psy.container)

class Config:
    def __init__(self):
        self.config = {}
        self.load_default_config()
    
    def load_default_config(self) -> None:
        # calc path
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_filename = os.path.join(script_dir, "config.jsonc")

        # load
        with open(config_filename, "r") as fp:
            config = json5.load(fp)

        # save
        self.config = config

    def select_file(self, file_in: str) -> bool:
        # get filename
        filename = os.path.basename(file_in)

         # check accept it
        if not filename in self.config['files']:
            return False

        # extract possible custom config
        return True

    def get(self, key: str, default = None):
        return self.config['default'].get(key, default)

def main(file_in: str, file_out:str) -> None:
    print(file_out)
    # vars
    config = Config()
    is_to_treat = config.select_file(file_in)

    # skip not to treat
    if not is_to_treat:
        shutil.copyfile(file_in, file_out)
        raise Exception(file_in)
        return

    # load
    reader = FortranReader()
    container = reader.psyir_from_file(file_in, free_form=False)

    # transform
    instumenter = TwinCheckerInstrumenter(file_in, file_out, config)
    instumenter.apply(container)

    # regen
    writer = FortranWriter()
    new_code = writer(container)

    # make line limites
    reshaper = FortLineLength()
    new_code_reshaped = reshaper.process(new_code)

    # replace __LINE__
    lines = new_code_reshaped.split("\n")
    for id, line in enumerate(lines):
        lines[id] = line.replace("__LINE__", str(id + 1))

    # reassemble
    code_final = "\n".join(lines)
    if not "twin_check" in code_final:
        raise Exception("lkj")

    # write down
    with open(file_out, "w+") as fp:
        fp.write(code_final)

if __name__ == "__main__":
    # common base options
    parser = argparse.ArgumentParser(description="Instumenter to inject twin_check calls in your fortran file.")

    # args
    parser.add_argument("-o", "--output", help="Output file to store result in.", required=True, type=str)
    parser.add_argument("input_file", help="The input file to parse and instrument.", type=str, nargs=1)

    # parse
    args = parser.parse_args()

    # call
    main(args.input_file[0], args.output)

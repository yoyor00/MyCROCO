#!/usr/bin/env python3

# python
import os
import json5
import shutil
import random
import argparse
import abc
# psyclone
from psyclone.psyir.nodes import Routine, Loop, Reference, Assignment, Call, BinaryOperation, UnaryOperation, Schedule, Literal, Node, IntrinsicCall, Container, IfBlock, Operation, ArrayReference, ACCDirective, RegionDirective
from psyclone.psyir.symbols import RoutineSymbol, NoType, ContainerSymbol, DataType, Symbol, ArrayType, CHARACTER_TYPE, INTEGER4_TYPE, INTEGER_UNDEF_TYPE, INTEGER8_TYPE, ScalarType
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.line_length import FortLineLength
from psyclone.f2pygen import DirectiveGen
# local
from trans_array_indices_to_get_full import ReferenceIndices2ArrayRangeTrans

##########################################################
class ACCCustomDirective(ACCDirective, RegionDirective, metaclass=abc.ABCMeta):
    '''
    Represent the "acc set devicenum()" directive.
    '''
    def __init__(self, parent=None, children=None, directive: str = ""):
        # A Directive always contains a Schedule
        super().__init__(children=children, parent=parent)
        self._directive = directive

    def __eq__(self, other):
        '''
        Checks whether two nodes are equal. Two ACCKernelsDirective nodes are
        equal if their default_present members are equal.

        :param object other: the object to check equality to.

        :returns: whether other is equal to self.
        :rtype: bool
        '''
        is_eq = super().__eq__(other)
        is_eq = is_eq and self._directive == other._directive

        return is_eq

    @property
    def directive(self):
        '''
        :returns: whether the "default(present)" clause is added to the \
                  kernels directive.
        :rtype: bool
        '''
        return self._directive

    def gen_code(self, parent):
        '''
        Generate the f2pygen AST entries in the Schedule for this
        OpenACC Kernels directive.

        :param parent: the parent Node in the Schedule to which to add this \
                       content.
        :type parent: sub-class of :py:class:`psyclone.f2pygen.BaseGen`

        '''
        self.validate_global_constraints()

        # We re-use the 'begin_string' method but must skip the leading 'acc'
        # that it includes.
        parent.add(DirectiveGen(parent, "acc", "begin",
                                *self.begin_string().split()[1:]))

        self.gen_post_region_code(parent)

    def begin_string(self):
        '''Returns the beginning statement of this directive, i.e.
        "acc kernels ...". The backend is responsible for adding the
        correct directive beginning (e.g. "!$").

        :returns: the beginning statement for this directive.
        :rtype: str

        '''
        result = self._directive

        return result

    def end_string(self):
        '''
        Returns the ending statement for this directive. The backend is
        responsible for adding the language-specific syntax that marks this
        as a directive.

        :returns: the closing statement for this directive.
        :rtype: str

        '''
        return ""


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
        self.enable_profile_after_kernel = self.config.get('profile_after_kernel', False)
        self.enable_profile_before_kernel = self.config.get('profile_before_kernel', False)
        self.enable_before_after_only_arrays = self.config.get('before_after_only_arrays', False)
        self.enable_profile_arrays = self.config.get('profile_arrays', False)
        self.enable_fixable = self.config.get('fixable', False)

        if ("ana_initial" in source_file):
            print("---------------------------")
            print(self.enable_profile_in_values)
            print(self.enable_profile_out_values)
            print(self.enable_profile_integers)
            print(self.enable_profile_after_kernel)
            print(self.enable_profile_before_kernel)
            print(self.enable_fixable)
            print("---------------------------")

    def gen_call(self, node: Node, value: Node = None) -> Call:
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
        if isinstance(node, Reference) and node.is_array and self.enable_fixable:
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
        if value != None:
            equation_str += " = " + value.debug_string()
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

    def gen_call_array(self, node: ArrayReference) -> Call:
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
        if isinstance(node, Reference) and node.is_array and self.enable_fixable:
            fixable = "_fixable"

        if "REAL" == t:
            twin_check_float = RoutineSymbol(f"twin_check_double{fixable}_array", NoType())
        elif "INTEGER" == t:
            twin_check_float = RoutineSymbol(f"twin_check_integer{fixable}_array", NoType())
        elif "BOOLEAN" == t:
            twin_check_float = RoutineSymbol(f"twin_check_bool{fixable}_array", NoType())
        else:
            raise Exception(str(literal))
            twin_check_float = RoutineSymbol("twin_check_float_array", NoType())
        equation_str = node.name + "(.....)"
        if 'zob' == equation_str:
            return None
        if ':' in equation_str:
            return None
        site_id = random.randint(0, 10000000000)
        self.site_ids.append(site_id)
        name = node.name
        reader = FortranReader()
        routine: Routine = node.ancestor(Routine)
        arguments = [
            reader.psyir_from_expression(f"{name}", routine.symbol_table),
            reader.psyir_from_expression(f"SIZE({name})", routine.symbol_table),
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
                loc_insert['call'] = self.gen_call(store_in, value)
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
        dest_file = self.dest_file #.replace('/home/svalat/Projects/minicroco/BENCH/rundir/BASIN-SMALL', './') #@todo: fix to keep full name but as isusue with line reshapingm
        routines = container.walk(Routine)
        enable_profile_after_kernel = self.enable_profile_after_kernel
        enable_profile_before_kernel = self.enable_profile_before_kernel

        for routine in routines:
            module_symbol = ContainerSymbol("twin_checker", True)
            routine.symbol_table.add(module_symbol)

            self.apply_on_node(routine, stop_type=Loop)
            for loop in routine.walk(Loop, stop_type=Loop):
                self.apply_on_node(loop)

            if enable_profile_before_kernel:
                for loop in routine.walk(Loop, stop_type=Loop):
                    parent = loop.parent
                    insert_pos = loop.position
                    acc_update_pos = insert_pos
                    set_seen = []
                    vars_to_update = []
                    for assign in loop.walk(Assignment):
                        set_seen.append(assign.lhs.debug_string())
                        # @todo make a single if many write to same var
                        for ref in assign.rhs.walk(Reference):
                            if "cff1_w" in set_seen:
                                print("-------------------")
                                print(ref.debug_string())
                                print(set_seen)
                            if ref.debug_string() not in set_seen and (not self.enable_before_after_only_arrays or isinstance(ref, ArrayReference)):
                                vars_to_update.append(ref.name)
                                if "cff1_w" in set_seen:
                                    print("took")
                                if isinstance(ref, ArrayReference):
                                    to_insert = {
                                        'parent': parent,
                                        'pos': insert_pos,
                                        'call': self.gen_call_array(ref)
                                    }
                                else:
                                    to_insert = {
                                        'parent': parent,
                                        'pos': insert_pos,
                                        'call': self.gen_call(ref)
                                    }
                                # insert
                                if to_insert['call'] != None:
                                    parent.addchild(to_insert['call'], to_insert['pos'])
                                    insert_pos += 1

                    # inject acc update
                    if len(vars_to_update) > 0 and False:
                        vars_to_update_str = ', '.join(set(vars_to_update))
                        dir_copy_to_cpu = ACCCustomDirective(parent=None, children=None, directive=f"acc update self ({vars_to_update_str})")
                        parent.addchild(dir_copy_to_cpu, acc_update_pos)
                        dir_copy_to_cpu = ACCCustomDirective(parent=None, children=None, directive=f"acc update device ({vars_to_update_str})")
                        parent.addchild(dir_copy_to_cpu, insert_pos)

            if enable_profile_after_kernel:
                for loop in routine.walk(Loop, stop_type=Loop):
                    parent = loop.parent
                    insert_pos = loop.position + 1
                    acc_update_pos = insert_pos
                    vars_to_update = []
                    for assign in loop.walk(Assignment):
                        # @todo make a single if many write to same var
                        ref: Reference = assign.lhs
                        if ref.name in vars_to_update:
                            continue
                        if self.enable_before_after_only_arrays and not isinstance(ref, ArrayReference):
                            continue
                        if isinstance(ref, ArrayReference) and not self.enable_profile_arrays:
                            continue
                        vars_to_update.append(ref.name)
                        if isinstance(ref, ArrayReference):
                            to_insert = {
                                'parent': parent,
                                'pos': insert_pos,
                                'call': self.gen_call_array(ref)
                            }
                        else:
                            to_insert = {
                                'parent': parent,
                                'pos': insert_pos,
                                'call': self.gen_call(ref)
                            }
                        # insert
                        if to_insert['call'] != None:
                            insert_pos += 1
                            parent.addchild(to_insert['call'], to_insert['pos'])

                    # inject acc update
                    if len(vars_to_update) > 0 and False:
                        vars_to_update_str = ', '.join(set(vars_to_update))
                        dir_copy_to_cpu = ACCCustomDirective(parent=None, children=None, directive=f"acc update self ({vars_to_update_str})")
                        parent.addchild(dir_copy_to_cpu, acc_update_pos)
                        dir_copy_to_cpu = ACCCustomDirective(parent=None, children=None, directive=f"acc update device ({vars_to_update_str})")
                        parent.addchild(dir_copy_to_cpu, insert_pos)

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

###########################################################
def simplify_file_name(filename: str) -> str:
    # calc cleaned file name (without .no-acc.cpp.mpc.......)
    # So keep [0].[-1]
    parts = filename.split('.')
    extention = parts[-1]
    name = parts[0]
    return f'{name}.{extention}'.replace('.F90', '.F')

class Config:
    def __init__(self):
        self.config = {}
        self.filename = 'UNDEFINED'
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
        simplified_name = simplify_file_name(filename)

         # check accept it
        if not simplified_name in self.config['files']:
            return False
        
        # keep track
        self.filename = filename

        # extract possible custom config
        return True

    def get(self, key: str, default = None):
        default_value = self.config['defaults'].get(key, default)
        for name, group in self.config.get('specific', {}).items():
            for file in group.get('files', []):
                if file == simplify_file_name(self.filename):
                    return group.get('config', {}).get(key, default_value)
        return default_value

def main(file_in: str, file_out:str) -> None:
    print(file_out)
    # vars
    config = Config()
    is_to_treat = config.select_file(file_in)

    # skip not to treat
    if not is_to_treat:
        raise Exception("File not supported !")
        shutil.copyfile(file_in, file_out)
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

    # split paths
    # in_string = False
    # is_path = False
    # splitted = ""
    # new_val = ""
    # for i, letter in enumerate(new_code):
    #     in_string = "'"
    #     escaped = False
    #     new_val = letter
    #     if letter == "\\":
    #         escaped = True
    #     elif letter == "/" and in_string and is_path:
    #         new_val = "/' '"
    #     elif letter == "'" and not escaped:
    #         in_string = not in_string
    #         if not in_string:
    #             is_path = False
    #         if new_code[i+1] == "/":
    #             is_path = True
    #     splitted += new_val
    # new_code = splitted

    # make line limites
    reshaper = FortLineLength()
    reshaper._key_lists['unknown'].append("/")
    reshaper._key_lists['statement'].append("/")
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

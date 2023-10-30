#!/usr/bin/env python3
##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import re
import os
import json

##########################################################
class CrocoVariable:
    def __init__(self, infos: dict):
        self.name = infos['name']
        self.type = infos['type']
        self.latex = infos['latex']
        self.longName = infos['longName']
        self.description = infos['description']
        self.docUrl = infos['docUrl']
        self.kind = infos['kind']
        self.usedAs = infos['usedAs']
        self.gridOrigin = infos['gridOrigin']
        self.compatVarName = infos['compatVarName']
        self.compatOrigFile = infos['compatOrigFile']
        self.sourceLine = infos['sourceLine']
        self.cppKeyCond = infos['cppKeyCond']

    def to_dict(self):
        return {
            'name': self.name,
            'type': self.type,
            'latex': self.latex,
            'longName': self.longName,
            'description': self.description,
            'docUrl': self.docUrl,
            'kind': self.kind,
            'usedAs': self.usedAs,
            'gridOrigin': self.gridOrigin,
            'compatVarName': self.compatVarName,
            'compatOrigFile': self.compatOrigFile,
            'cppKeyCond': self.cppKeyCond,
            'sourceLine': self.sourceLine
        }

##########################################################
class CrocoVariableExtractor:
    def __init__(self):
        self.variables = {}

    def parse_file(self, fname: str) -> None:
        # parse macros for vars
        pattern_pos = r"^\s*!POSEIDON\s+([a-zA-Z0-9]+)\s*:\s*([a-zA-Z0-9():#\\/,._ -]+)(\s*!.*)?$"
        pattern_decl = r"^\s*(real|integer)\s+([a-zA-Z0-9_]+)\s*(.*)$"

        # extract condition of activation
        pattern_macro_if = r"^#\s*if (.*)\s*$"
        pattern_macro_ifdef = r"^#\s*ifdef (.*)\s*$"
        pattern_macro_else = r"^#\s*else.*$"
        pattern_macro_endif = r"^#\s*endif.*$"

        # open file
        with open(fname, 'r') as fp:
            # load all lines
            lines = fp.readlines()

            # state
            state_infos = {}
            state_macro_if_stack = []

            # loop on lines
            for line_id, line in enumerate(lines):
                # apply regexps
                capture_poseidon = re.match(pattern_pos, line)
                capture_decl = re.match(pattern_decl, line)

                # apply regexps
                capture_macro_if = re.match(pattern_macro_if, line)
                capture_macro_ifdef = re.match(pattern_macro_ifdef, line)
                capture_macro_else = re.match(pattern_macro_else, line)
                capture_macro_endif = re.match(pattern_macro_endif, line)

                # check for line error
                if "!POSEIDON" in line and not capture_poseidon:
                    raise Exception(f"Invalid poseidon line : '{line}'")
                
                # handle if/else/ifdef/endif
                if capture_macro_if:
                    state_macro_if_stack.append(capture_macro_if.group(1))
                    print(state_macro_if_stack)
                if capture_macro_ifdef:
                    state_macro_if_stack.append("defined "+capture_macro_ifdef.group(1))
                    print(state_macro_if_stack)
                if capture_macro_else:
                    state_macro_if_stack[-1] = "! ( " + state_macro_if_stack[-1] + " ) "
                    print(state_macro_if_stack)
                if capture_macro_endif:
                    state_macro_if_stack.pop()
                    print(state_macro_if_stack)

                # capture poseidon key/value
                if capture_poseidon:
                    print(capture_poseidon.group())

                    # extract
                    key = capture_poseidon.group(1)
                    value= capture_poseidon.group(2)

                    # display
                    print(f"{key}: {value}")

                    # store key/value for latter use (when will see the decl)
                    state_infos[key] = value

                # capture var decl (so we can flush)
                if capture_decl and len(state_infos.keys()) > 0:
                    print(capture_decl.group())

                    # extract
                    vtype = capture_decl.group(1)
                    vname = capture_decl.group(2)

                    # some coherency checks
                    assert vtype == state_infos['type']
                    assert vname == state_infos['compatVarName']
                    assert os.path.basename(fname) == state_infos['compatOrigFile']

                    # store #if stack
                    print(state_macro_if_stack)
                    state_infos["cppKeyCond"] = state_macro_if_stack.copy()

                    # line
                    state_infos["sourceLine"] = line_id

                    # push variable to the list
                    var = CrocoVariable(state_infos)
                    self.variables[var.name] = var
                    state_infos = {}

            # check that we were rightly capturing the #if / #ifdef / #else / #endif so it closes well
            assert len(state_macro_if_stack) == 0

    def to_dict(self):
        result = {}
        for key, var in self.variables.items():
            result[var.name] = var.to_dict()
        return result

##########################################################
if __name__ == "__main__":
    extractor = CrocoVariableExtractor()
    extractor.parse_file("../OCEAN/ocean3d.h")
    print(json.dumps(extractor.to_dict(), indent='\t'))
    with open("croco_variables.json", "w+") as fp:
        json.dump(extractor.to_dict(), indent='\t', fp=fp)
    with open("croco_variables.js", "w+") as fp:
        fp.write("var crocoVariables = ")
        json.dump(extractor.to_dict(), indent='\t', fp=fp)
        fp.write(";")

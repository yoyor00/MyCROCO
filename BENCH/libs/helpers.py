##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''Provide some helper function for example to run commands and patch files.'''

##########################################################
import re
import os
import sys
import shutil
import subprocess
import tempfile
import traceback
from termcolor import cprint
from timeit import timeit
from .messaging import Messaging
from contextlib import contextmanager

##########################################################
def display_run_error(command: str, outpout: str) -> None:
    dir = os.getcwd()
    keep_n_last_lines = 32
    truncated_output = '\n'.join(outpout.split('\n')[-keep_n_last_lines:])
    cprint("-----------------------------------------------", 'red')
    cprint(truncated_output, 'dark_grey')
    cprint("-----------------------------------------------", 'red')
    cprint(f"Error from command : {command}", 'red')
    cprint(f"Error from workdir : {dir}", 'red')
    cprint("-----------------------------------------------", 'red')
    print("")

##########################################################
def print_exception(exception: Exception) -> None:
    dir = os.getcwd()
    cprint("-----------------------------------------------", 'red')
    cprint(''.join(traceback.format_exception(exception)), 'dark_grey')
    cprint("-----------------------------------------------", 'red')
    cprint(f"Error from command : {' '.join(sys.argv)}", 'red')
    cprint(f"Error from workdir : {dir}", 'red')
    cprint("-----------------------------------------------", 'red')
    cprint(str(exception), 'red')
    cprint("-----------------------------------------------", 'red')
    print("")

##########################################################
def run_shell_command(command, capture = True, show_on_error = True):
    '''Print the command and run it. On failure it prints the output.'''

    # print command
    Messaging.command(command)

    # capture or not capture
    if capture:
        # run
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # on failure
        if result.returncode != 0:
            if show_on_error:
                display_run_error(command, result.stdout)
            raise Exception("Fail to run command !")
    else:
        print("-----------------------------------------------")
        subprocess.run(command, shell=True, check=True)
        print("-----------------------------------------------")

def run_shell_command_time(command, verbose: bool = False):
    '''Print the command and run it by measuring time. On failure it prints the output.'''
    # capture or not capture
    if verbose:
        return timeit(stmt = f"subprocess.run('{command}', check=True, shell=True)", setup = "import subprocess", number = 1)
    else:
        with tempfile.NamedTemporaryFile("w+") as log_fp:
            try:
                return timeit(stmt = f"subprocess.run('{command} 2>&1 > {log_fp.name}', check=True, shell=True)", setup = "import subprocess", number = 1)
            except:
                display_run_error(command, log_fp.read())
                raise Exception("Fail to run command !")

def replace_in_file(in_path: str, out_path: str, pattern: str, replace_by:str):
    '''Patch a file by replacing some patterns inside.'''

    # print progress
    fname = os.path.basename(out_path)
    Messaging.step(f"Patching {fname}")

    # open & load
    with open(in_path, 'r') as fp:
        content = fp.read()

    # patch
    content = content.replace(pattern, replace_by)

    # rewrite
    with open(out_path, 'w+') as fp:
        fp.write(content)

def contains_one_of(seach_in:str, elements) -> bool:
    '''check is the string contains one of the given elements.'''

    # if search a simple string
    if isinstance(elements, str):
        return elements in seach_in
    elif isinstance(elements, list):
        for element in elements:
            if element in seach_in:
                return True
        return False
    else:
        raise Exception("Unsupported type !")

def patch_lines(path: str, rules: list, allow_already_done = False):
    '''
    Patch a file by inserting some lines at the given place. The changes are
    requested via a list of rules to apply.
    '''
    fname = os.path.basename(path)

    # copy to keep old version
    shutil.copy(path, f'{path}.backup')

    # load
    with open(path, 'r') as fp:
        lines = fp.readlines()

    # loop on rules
    for rule in rules:
        # Display
        if 'descr' in rule:
            descr = rule['descr']
            Messaging.step(f"Patching {fname} [ {descr} ]")
        else:
            Messaging.step(f"Patching {fname}")
        # replace rule
        if rule['mode'] == 'replace':
            if not (allow_already_done and rule['by'] in lines):
                id = lines.index(rule['what'])
                lines[id] = rule['by']
        # insert after a given line
        elif rule['mode'] == 'insert-after':
            id = lines.index(rule['what'])
            if isinstance(rule['insert'], str):
                lines.insert(id + 1, rule['insert'])
            else:
                for i, line in enumerate(rule['insert']):
                    lines.insert(id + 1 + i, line)
        # insert before a given line
        elif rule['mode'] == 'insert-before':
            id = lines.index(rule['what'])
            if isinstance(rule['insert'], str):
                lines.insert(id, rule['insert'])
            else:
                for i, line in enumerate(rule['insert']):
                    lines.insert(id + i, line)
        # search a start line, go until a stop and insert before.
        elif rule['mode'] == 'insert-after-before':
            id = 0
            while not rule['after'] in lines[id]:
                id += 1
            id += 1
            while not contains_one_of(lines[id], rule['before']):
                id += 1
            lines.insert(id, rule['insert'])
        elif rule['mode'] == 'insert-at-end':
            lines.append(rule['what'])

    # save
    with open(path, 'w+') as fp:
        fp.writelines(lines)

@contextmanager
def move_in_dir(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

def apply_vars_in_str(value: str, vars: {}) -> str:
    # TODO: do not yet allow recursion but can ben done easly the way it is.

    # extract {case.value} patterns
    patterns = re.findall(r"\{[a-zA-Z0-9._-]+\}", value)

    # replace all
    if patterns:
        for pattern in patterns:
            # progress
            Messaging.step(f"Apply variable {pattern}...")

            # fetch value
            path = pattern.replace("{", "").replace('}','').split('.')

            # start
            cursor = vars
            for element in path:
                cursor = cursor[element]

            # replace
            value = value.replace(pattern, cursor)

    # ok
    return value

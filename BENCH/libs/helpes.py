##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''Provide some helper function for example to run commands and patch files.'''

##########################################################
import re
import os
import shutil
import subprocess
from timeit import timeit
from .messaging import Messaging
from contextlib import contextmanager

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
                print("-----------------------------------------------")
                print(result.stdout)
                print("-----------------------------------------------")
                print(f"Error from command : {command}")
                print("-----------------------------------------------")
            raise Exception("Fail to run command !")
    else:
        subprocess.run(command, shell=True, check=True)

def run_shell_command_time(command, verbose: bool = False):
    '''Print the command and run it by measuring time. On failure it prints the output.'''
    # capture or not capture
    if verbose:
        return timeit(stmt = f"subprocess.run('{command}', check=True, shell=True)", setup = "import subprocess", number = 1)
    else:
        return timeit(stmt = f"subprocess.run('{command} > /dev/null', check=True, shell=True)", setup = "import subprocess", number = 1)

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

def patch_lines(path: str, rules: list):
    '''
    Patch a file by inserting some lines at the given place. The changes are
    requested via a list of rules to apply.
    '''
    fname = os.path.basename(path)
    Messaging.step(f"Patching {fname}")

    # copy to keep old version
    shutil.copy(path, f'{path}.backup')

    # load
    with open(path, 'r') as fp:
        lines = fp.readlines()

    # loop on rules
    for rule in rules:
        # replace rule
        if rule['mode'] == 'replace':
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
    # extract {case.value} patterns
    patterns = re.findall(r"\{[a-zA-Z0-9._-]+\}", value)

    # replace all
    for pattern in patterns:
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

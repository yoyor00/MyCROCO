#!/usr/bin/env python3
##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
This script simply aims (for now) on removing the ACC
from some files

TODO
----
This can be removed when work will be finished
'''

###########################################################
import os
import json
import shutil
import argparse
from common.wrapper_helpers import get_rules_for_file

###########################################################
def parse_command_line():
    '''
    Parse the command line and return the user options.
    '''

    # build arg parsor
    parser = argparse.ArgumentParser(
        prog='psyclone.preprocessor.skip.acc.py',
        description='Patch the files to make them skiping OpenACC'
    )

    # define args
    parser.add_argument('IN_FILE', help="Input source file")
    parser.add_argument('OUT_FILE', help="Output file")

    # parser
    options = parser.parse_args()

    # ok
    return options

###########################################################
def patch_cppdef(in_file: str, out_file: str):
    '''
    Patch the cppdef.h file to disable OPENACC
    '''
    with open(in_file, 'r') as fp_in:
        # load content
        data = fp_in.read()

        # replace the definition of OPENACC keys in cppdef
        data = data.replace('# define OPENACC', '# undef OPENACC')

        # for the CMake version, we need to override and erase the definition done in  config_post.h
        data = data.replace('#include \"config_post.h\"', '#include \"config_post.h\"\n# undef OPENACC')

        # write again
        with open(out_file, 'w+') as fp_out:
            fp_out.write(data)

###########################################################
def patch_source_file(in_file: str, out_file: str, new_cppdef: str):
    '''
    Patch the source file to use the customized cppdef.h without OPENACC
    '''

    with open(in_file, 'r') as fp_in:
        # load content
        data = fp_in.read()

        # replace
        data = data.replace('cppdefs.h', new_cppdef)

        # write again
        with open(out_file, 'w+') as fp_out:
            fp_out.write(data)

###########################################################
if __name__ == '__main__':
    # get options
    options = parse_command_line()

    # extract what is usefulll
    in_file = options.IN_FILE
    out_file = options.OUT_FILE

    # compute some names
    in_file_name = os.path.basename(in_file)
    in_file_dir = os.path.dirname(in_file)

    # load config
    skip_acc = get_rules_for_file(in_file_name)['skip-acc']

    # if need to skip ACC
    if skip_acc:
        # patch cppdef and make this one for the file
        cppdef = os.path.join(in_file_dir, 'cppdefs.h')
        cppdef_patched = f"{out_file}.cppdefs-no-acc.h"
        patch_cppdef(cppdef, cppdef_patched)

        # patch the source file to move if to its own cppdef without OPENACC
        patch_source_file(in_file, out_file, cppdef_patched)
    else:
        shutil.copy(in_file, out_file)

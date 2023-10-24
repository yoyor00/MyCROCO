#!/usr/bin/env python3
##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

import os
import sys
import json
import subprocess
from common.wrapper_helpers import get_rules_for_file, get_local_file_path, extract_source_file_from_args

###########################################################
'''
This script aims at wrapping the compiler to apply psyclone
with the given transformation on the files listed in the
config file.

The way to use is is just giving the same options you would
have given to the compiler and just define the real compiler
to be used by adding the compiler as first option.

USAGE :
  ./psyclone.compiler.wrapper.py gcc -I./src my_file.F90 -o my_file.o
'''

###########################################################
if __name__ == '__main__':
    # default
    transformation_script = None

    # check if compiling a .o file
    is_for_o_file = ('-c' in sys.argv)

    # if possibly needs psyclone
    if is_for_o_file:
        # extract source file
        source_file = extract_source_file_from_args(sys.argv)
        assert source_file is not None

        # get path of current script
        current_script = os.path.realpath(__file__)
        script_path = os.path.dirname(current_script)

        # Get transformation script
        transformation_script = get_rules_for_file(source_file)['script']

        # call psyclone if needed
        if transformation_script is not None:
            # Generate the dummy version to be able to diff easily to see what we injected
            # Remark:  This is only for debugging purpose
            psyclone_dummy_source_file = source_file.replace(".F", ".psyclone.dummy.F90")
            subprocess.run([
                            'psyclone',
                            '-api', 'nemo',
                            '-l' , 'output',
                            '-s', get_local_file_path(os.path.join('scripts', 'trans_dummy.py')),
                            '-opsy' ,psyclone_dummy_source_file, source_file
                        ],
                        check=True)

            # run psyclone
            psyclone_source_file = source_file.replace(".F", ".psyclone.F90")
            subprocess.run([
                            'psyclone',
                            '-api', 'nemo',
                            '-l' , 'output',
                            '-s', get_local_file_path(os.path.join('scripts', transformation_script)),
                            '-opsy' ,psyclone_source_file, source_file
                        ],
                        check=True)

            # replace source filestart
            pos = sys.argv.index(source_file)
            sys.argv[pos] = psyclone_source_file

    # call compiler
    cmd = sys.argv[1:]
    subprocess.run(cmd, check=True)

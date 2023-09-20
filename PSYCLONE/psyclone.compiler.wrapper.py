#!/usr/bin/env python3
##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

import os
import sys
import subprocess

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
def extract_source_file(args: list):
    '''
    Extract the fortran source file from the compiler command line.
    
    Remark: we consider here a single fortran file compile at a time.

    Parameters
    ----------
    args: list
        The argument list in which to search.

    Returns
    -------
    Path of the source file we are compiling.
    '''

    # search first
    for value in args:
        if value.endswith(('.F', '.f', '.F90', '.f90')):
            return value
    
    # not found
    raise Exception("No fortran file found in command line !")

###########################################################
def get_transformation_for_file(config_file: str, source_file: str):
    '''
    Extract the transformation script to apply to the given file
    from the textual config file given as parameter.

    Parameters
    ----------
    config_file : str
        Path of the textual config file to load.
    source_file : str
        The current source file to compile and for which we want
        to extract the psyclone script to use.

    Returns
    -------
    Path to the script to be used or None.
    '''

    # extract file name
    file_name = os.path.basename(source_file)

    # calc cleaned file name (without .no-acc.cpp.mpc.......)
    file_name_simple = '.'.join([file_name.split('.')[i] for i in [0,-1]])

    # search in config file
    with open(config_file, 'r') as fp:
        # load
        lines = fp.readlines()
        for line in lines:
            entries = line.replace('\n','').split('\t')
            if entries[0] == file_name or entries[0] == file_name_simple:
                return entries[1]

    # not found
    return None

###########################################################
if __name__ == '__main__':
    # default
    transformation_script = None

    # check if compiling a .o file
    is_for_o_file = ('-c' in sys.argv)

    # if possibly needs psyclone
    if is_for_o_file:
        # extract source file
        source_file = extract_source_file(sys.argv)
        assert source_file is not None

        # get path of current script
        current_script = os.path.realpath(__file__)
        script_path = os.path.dirname(current_script)

        # Get transformation script
        config_file = os.path.join(script_path, "psyclone.rules.lst")
        transformation_script = get_transformation_for_file(config_file, source_file)

        # call psyclone if needed
        if transformation_script is not None:
            # Generate the dummy version to be able to diff easily to see what we injected
            # Remark:  This is only for debugging purpose
            psyclone_dummy_source_file = source_file.replace(".F", ".psyclone.dummy.F90")
            subprocess.run([
                            'psyclone',
                            '-api', 'nemo',
                            '-l' , 'output',
                            '-s', os.path.join(script_path, 'scripts', 'trans_dummy.py'),
                            '-opsy' ,psyclone_dummy_source_file, source_file
                        ],
                        check=True)

            # run psyclone
            psyclone_source_file = source_file.replace(".F", ".psyclone.F90")
            subprocess.run([
                                'psyclone',
                                '-api', 'nemo',
                                '-l' , 'output',
                                '-s', os.path.join(script_path, 'scripts', transformation_script),
                                '-opsy' ,psyclone_source_file, source_file
                            ],
                            check=True)

            # replace source filestart
            pos = sys.argv.index(source_file)
            sys.argv[pos] = psyclone_source_file

    # call compiler
    cmd = sys.argv[1:]
    subprocess.run(cmd, check=True)

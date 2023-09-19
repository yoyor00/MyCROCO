#!/usr/bin/env python3

import os
import sys
import subprocess

###########################################################
# This script aims at wrapping the compiler to apply psyclone
# with the given transformation on the files listed in the
# config file.
#
# The way to use is is just giving the same options you would
# have given to the compiler and just define the real compiler
# to be used by adding the compiler as first option.

###########################################################
def extract_compiler_and_remove_from_list(args: list):
    '''Search for compiler name by searching --compiler option'''
    # search and extract
    compiler = args[1]

    # remove from options to get only what to forward to compiler
    del args[1]

    # ok
    return compiler

###########################################################
def extract_source_file(args: list):
    source = ''
    for value in args:
        if value.endswith(('.F', '.f', '.F90', '.f90')):
            source = value
    return source

###########################################################
def get_transformation_for_file(config_file: str, source_file: str):
    # extract file name
    file_name = os.path.basename(source_file)

    # search in config file
    with open(config_file, 'r') as fp:
        # load
        lines = fp.readlines()
        for line in lines:
            entries = line.replace('\n','').split('\t')
            if entries[0] == file_name:
                return entries[1]
    
    # not found
    return None

###########################################################
if __name__ == '__main__':
    # extracct compiler from command line
    real_compiler = extract_compiler_and_remove_from_list(sys.argv)

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
    cmd = [real_compiler] + sys.argv[1:]
    #print(' '.join(cmd))
    subprocess.run(cmd, check=True)
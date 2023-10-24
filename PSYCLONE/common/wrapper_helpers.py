##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json

###########################################################
def get_local_file_path(file_name: str):
    '''
    Calculate the absolute path of file in the PSYCLONE directory
    '''

    # find PSYCLONE dir path
    current_script = os.path.realpath(__file__)
    psyclone_dir_path = os.path.dirname(os.path.dirname(current_script))

    # Get transformation script
    path = os.path.join(psyclone_dir_path, file_name)

    # ok
    return path

###########################################################
def simplify_file_name(filename: str) -> str:
    # calc cleaned file name (without .no-acc.cpp.mpc.......)
    # So keep [0].[-1]
    parts = filename.split('.')
    extention = parts[-1]
    name = parts[0]
    return f'{name}.{extention}'

###########################################################
def get_rules_for_file(source_file: str) -> dict:
    '''
    Extract the transformation script to apply to the given file
    from the textual config file given as parameter.
    '''

    # get config file path
    config_file = get_local_file_path('psyclone.rules.json')

    # extract file name
    file_name = os.path.basename(source_file)

    # calc cleaned file name (without .no-acc.cpp.mpc.......)
    file_name_simple = simplify_file_name(file_name)

    # default
    rule = {
        'script': None,
        'skip-acc': False
    }

    # open config file
    with open(config_file, 'r') as fp:
        # load as JSON
        config = json.load(fp)

        # perform some checks
        assert 'rules' in config
        rules = config['rules']
        assert isinstance(config['rules'], dict)

        # load
        if file_name_simple in rules:
            # extract
            file_rule = rules[file_name_simple]

            # apply default when needed
            for key, value in file_rule.items():
                rule[key] = value

    # debug
    #print({file_name_simple: rule})

    # not found
    return rule

###########################################################
def extract_source_file_from_args(args: list):
    '''
    Extract the fortran source file from the compiler command line.
    
    Remark: we consider here a single fortran file compile at a time.
    '''

    # search first
    for value in args:
        if value.endswith(('.F', '.f', '.F90', '.f90')):
            return value
    
    # not found
    raise Exception("No fortran file found in command line !")

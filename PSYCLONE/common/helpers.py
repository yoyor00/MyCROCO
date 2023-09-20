##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
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
def get_rules_for_file(source_file: str):
    '''
    Extract the transformation script to apply to the given file
    from the textual config file given as parameter.
    '''

    # get config file path
    config_file = get_local_file_path('psyclone.rules.json')

    # extract file name
    file_name = os.path.basename(source_file)

    # calc cleaned file name (without .no-acc.cpp.mpc.......)
    # So keep [0].[-1]
    parts = file_name.split('.')
    extention = parts[-1]
    name = parts[0]
    file_name_simple = f'{name}.{extention}'

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
            for key, value in enumerate(file_rule):
                rule[key] = value

    # not found
    return rule

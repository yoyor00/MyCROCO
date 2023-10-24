##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some unit tests about the wrapper helper functions.
'''

##########################################################
# internal
from common.wrapper_helpers import simplify_file_name, get_rules_for_file, extract_source_file_from_args

##########################################################
def test_simplify_file_name():
    assert simplify_file_name('simple.F') == 'simple.F'
    assert simplify_file_name('simple.no-acc.mpc.F') == 'simple.F'

##########################################################
def test_get_rules_for_file():
    assert get_rules_for_file('fake_for_unit_test.acc.F90') == {"script": "fake_script.py", "skip-acc": True}

##########################################################
def test_extract_source_file_from_args():
    assert extract_source_file_from_args(['-Wall', '-I./tmp', '-o','tmp.o', 'tmp.F90']) == 'tmp.F90'
    assert extract_source_file_from_args(['-Wall', '-I./tmp', '-o','tmp.o', 'tmp.F']) == 'tmp.F'
    assert extract_source_file_from_args(['-Wall', '-I./tmp', '-o','tmp.o', 'tmp.f90']) == 'tmp.f90'
    assert extract_source_file_from_args(['-Wall', '-I./tmp', '-o','tmp.o', 'tmp.f']) == 'tmp.f'

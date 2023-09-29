##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation helper functins.
'''

##########################################################
from helper import *
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader

##########################################################
def test_add_missing_device_vars_step3d():
    '''
    Check whether it inserts the vars in step3d_t.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(Istr,Iend,Jstr,Jend)
            implicit none
            integer*4 Istr,Iend,Jstr,Jend, i,j
            integer*4 IstrR,IendR,JstrR,JendR
            if (istr.eq.1) then
                IstrR=Istr-1
            endif
        end
    ''', free_form = True)

    # apply
    add_missing_device_vars(root_node)

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert 'integer, parameter :: my_acc_device = 0' in gen_source
    assert 'logical, parameter :: compute_on_device = .true.' in gen_source

##########################################################
def test_add_missing_device_vars_not_apply():
    '''
    Check whether it not inserts in others.
    '''

    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine dummy(Istr,Iend,Jstr,Jend)
            implicit none
            integer*4 Istr,Iend,Jstr,Jend, i,j
            integer*4 IstrR,IendR,JstrR,JendR
            if (istr.eq.1) then
                IstrR=Istr-1
            endif
        end
    ''', free_form = True)

    # apply
    add_missing_device_vars(root_node)

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    assert 'my_acc_device' not in gen_source
    assert 'compute_on_device' not in gen_source

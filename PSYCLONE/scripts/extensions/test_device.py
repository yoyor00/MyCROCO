##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper functions.
'''

##########################################################
from .device import *
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
        subroutine step3d_t(istr)
            implicit none
            integer*4 istr
            integer*4 istrr
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
    assert gen_source == '''\
subroutine step3d_t(istr)
  integer, parameter :: my_acc_device = 0
  logical, parameter :: compute_on_device = .true.
  integer*4 :: istr
  integer*4 :: istrr

  if (istr == 1) then
    istrr = istr - 1
  end if

end subroutine step3d_t
'''

##########################################################
def test_add_missing_device_vars_not_apply():
    '''
    Check whether it not inserts in others.
    '''

    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine dummy(istr)
            implicit none
            integer*4 istr
            integer*4 istrr
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
    assert gen_source == '''\
subroutine dummy(istr)
  integer*4 :: istr
  integer*4 :: istrr

  if (istr == 1) then
    istrr = istr - 1
  end if

end subroutine dummy
'''

##########################################################
def test_set_device_type():
    '''
    Check if well inject the expected code.
    '''

    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(istr)
            implicit none
            integer*4 istr
            integer*4 tile
            call child_subroutine(istr)
        end
    ''', free_form = True)

    # apply
    set_device_tile(root_node)

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert gen_source == '''\
subroutine step3d_t(istr)
  integer*4 :: istr
  integer*4 :: tile

  !$acc set device_num(tile)
  call child_subroutine(istr)

end subroutine step3d_t
'''

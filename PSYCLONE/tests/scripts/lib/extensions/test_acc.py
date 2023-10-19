##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper functions.
'''

##########################################################
# psyclone
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.transformations import ACCParallelTrans
# internal
from scripts.lib.extensions.acc import *

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

##########################################################
def test_set_private_on_loop():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real tt(N, 10, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        tt(i, j, k) = i+j+k
                    end do
                end do
                do j = 1, 10
                    do i = 1, 20
                        tt(i, j, k) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    k_loop = root_node.walk(Loop)[0]
    assert k_loop.variable.name == 'k'

    # patch
    ACCParallelTrans().apply(k_loop, options={})
    set_private_on_loop(k_loop, 'i', ['tt'])

    # regen
    gen_source = FortranWriter()(root_node)
    print(gen_source)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  real, dimension(n,10,20) :: tt

  !$acc parallel default(present)
  do k = 1, 100, 1
    do j = 1, 10, 1
      !$acc loop independent private(tt)
      do i = 1, 20, 1
        tt(i,j,k) = i + j + k
      enddo
    enddo
    do j = 1, 10, 1
      !$acc loop independent private(tt)
      do i = 1, 20, 1
        tt(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

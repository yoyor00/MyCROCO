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
# internal
from scripts.extensions.scratch import *

##########################################################
def test_add_1d_scratch_var():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N
            implicit none
        end
    ''', free_form = True)

    # apply
    routine = root_node.walk(Routine)[0]
    add_1d_scratch_var(routine, 'dc')

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  real, dimension(0:n) :: dc1d


end subroutine step3d_t
'''

##########################################################
def test_patch_scratch_1d_arrays():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j
            real dc(N, 10)
            implicit none
            do i = 1, 10
                do j = 1, 10
                    dc(i, j) = i+j
                end do
            end do
        end
    ''', free_form = True)

    # apply
    routine = root_node.walk(Routine)[0]
    add_1d_scratch_var(routine, "dc")
    patch_scratch_1d_arrays(routine, ['dc'])

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  real, dimension(n,10) :: dc
  real, dimension(0:n) :: dc1d

  do i = 1, 10, 1
    do j = 1, 10, 1
      dc1d(j) = i + j
    enddo
  enddo

end subroutine step3d_t
'''

##########################################################
def test_add_3d_scratch_var():
    '''
    Check whether it inserts the vars subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t()
            integer*4 N
            integer*4 A3d(10,10,10)
            call step3d_tile(N)
        end
        subroutine step3d_tile(N)
            integer*4 N
            integer*4 istr, iend, jstr, jend
            implicit none
        end
    ''', free_form = True)

    # apply
    routine = root_node.walk(Routine)[1]
    assert routine.name == 'step3d_tile'
    add_3d_scratch_var(root_node, routine, 'dc', 1)

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert gen_source == '''\
subroutine step3d_t()
  integer*4 :: n
  integer*4, dimension(10,10,10) :: a3d

  call step3d_tile(n, a3d(1,1,trd))

end subroutine step3d_t
subroutine step3d_tile(n, dc_3d)
  integer*4 :: n
  real, dimension(istr - 2:iend + 2,jstr - 2:jend + 2,n), intent(inout) :: dc_3d
  integer*4 :: istr
  integer*4 :: iend
  integer*4 :: jstr
  integer*4 :: jend


end subroutine step3d_tile
'''

##########################################################
def test_patch_scratch_3d_arrays():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            integer*4 istr, jstr, iend, jend
            real fx(N, 20)
            implicit none
            do i = 1, 10
                do j = 1, 20
                    do k = 1, 10
                        fx(i, j) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    routine = root_node.walk(Routine)[0]
    add_3d_scratch_var(root_node, routine, "fx", 1)
    patch_scratch_3d_arrays(routine, ['fx'])

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    print(gen_source)
    assert gen_source == '''\
subroutine step3d_t(n, fx_3d)
  integer*4 :: n
  real, dimension(istr - 2:iend + 2,jstr - 2:jend + 2,n), intent(inout) :: fx_3d
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  integer*4 :: istr
  integer*4 :: jstr
  integer*4 :: iend
  integer*4 :: jend
  real, dimension(n,20) :: fx

  do i = 1, 10, 1
    do j = 1, 20, 1
      do k = 1, 10, 1
        fx_3d(i,j,k) = i + j + k
      enddo
    enddo
  enddo

end subroutine step3d_t
'''

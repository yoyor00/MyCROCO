##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper
functions.
'''

##########################################################
import pytest
from ..loops_helpers import *
from psyclone.psyir.nodes import Node
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader

##########################################################
def test_extract_loop_indices_order():
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
    first_loop = root_node.walk(Loop)[0]
    order = extract_loop_indices_order(first_loop)
    
    # check
    assert order == ['i', 'j', 'k']

##########################################################
def test_get_first_loop_on():
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
                ! first J loop to capture
                do j = 1, 20
                    ! first K loop to capture
                    do k = 1, 10
                        fx(i, j) = i+j+k
                    end do
                end do
                do j = 1, 20
                    fx(i, j) = i+j
                end do
            end do
        end
    ''', free_form = True)

    # apply
    first_loop = root_node.walk(Loop)[0]
    k_loop = get_first_loop_on(first_loop, "k")
    first_j_loop = get_first_loop_on(first_loop, "j")
    
    # check
    assert k_loop.variable.name == 'k'
    assert first_j_loop.variable.name == 'j'

    # check if got the first one, which contains a second loop.
    assert len(first_j_loop.walk(Loop)) == 2

##########################################################
def test_detach_and_get_childs():
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
                ! the content we want to detach
                fx(i, j) = i+j
                fx(i+1, j) = 10
                fx(i-1, j) = -10
            end do
            do k = 1, 20
                ! where we want to inject it
            end do
        end
    ''', free_form = True)

    # get the first loop
    i_loop = root_node.walk(Loop)[0]

    # extract its content
    i_loop_detatched_childs = detach_and_get_childs(i_loop)

    # re-inject it inside the other loop
    k_loop = root_node.walk(Loop)[1]
    for child in i_loop_detatched_childs:
        k_loop.loop_body.children.append(child)

    # remove i loop (which is empty)
    assert len(i_loop.loop_body.children) == 0
    #i_loop.detach()

    # regen
    gen_source = FortranWriter()(root_node)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  integer*4 :: istr
  integer*4 :: jstr
  integer*4 :: iend
  integer*4 :: jend
  real, dimension(n,20) :: fx

  do i = 1, 10, 1
  enddo
  do k = 1, 20, 1
    fx(i,j) = i + j
    fx(i + 1,j) = 10
    fx(i - 1,j) = -10
  enddo

end subroutine step3d_t
'''

##########################################################
def is_loop_using_var():
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
    first_loop = root_node.walk(Loop)[0]
    k_loop = get_first_loop_on(first_loop, "k")

    # check
    assert is_loop_using_var(k_loop, "k")
    assert not is_loop_using_var(k_loop, "l")

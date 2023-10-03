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
from ..loops import *
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.transformations import ACCParallelTrans

##########################################################
def test_handle_kji_loop_single():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    k_loop = root_node.walk(Loop)[0]
    assert k_loop.variable.name == 'k'

    # patch
    handle_kji_loop(k_loop, ['fx'])

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
  real, dimension(n,20) :: fx

  do k = 1, 100, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        fx_3d(i,j,k) = i + j + k
      enddo
    enddo
  enddo

end subroutine step3d_t
'''

##########################################################
def test_handle_kji_loop_splitting():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
                do j = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    k_loop = root_node.walk(Loop)[0]
    assert k_loop.variable.name == 'k'

    # patch
    handle_kji_loop(k_loop, ['fx'])

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
  real, dimension(n,20) :: fx

  do k = 1, 100, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        fx_3d(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  do k = 1, 100, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        fx_3d(i,j,k) = i + j + k
      enddo
    enddo
  enddo

end subroutine step3d_t
'''

##########################################################
def test_handle_kji_loop_splitting_with_if():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
                if (N == 55) then
                    do j = 1, 10
                        do i = 1, 20
                            fx(i, j) = i+j+k
                        end do
                    end do
                end if
            end do
        end
    ''', free_form = True)

    # apply
    k_loop = root_node.walk(Loop)[0]
    assert k_loop.variable.name == 'k'

    # patch
    handle_kji_loop(k_loop, ['fx'])

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
  real, dimension(n,20) :: fx

  do k = 1, 100, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        fx_3d(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  if (n == 55) then
    do k = 1, 100, 1
      do j = 1, 10, 1
        do i = 1, 20, 1
          fx_3d(i,j,k) = i + j + k
        enddo
      enddo
    enddo
  end if

end subroutine step3d_t
'''

##########################################################
def test_handle_kji_loop_splitting_with_if_conflict():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
                !---------------------------------------------------------------
                ! by default the function currently try to make the k-loop
                ! crossing this line which is invalid.
                ! Not seen in current CROCO, but might handle if it happens
                !---------------------------------------------------------------
                if (k == 55) then
                    do j = 1, 10
                        do i = 1, 20
                            fx(i, j) = i+j+k
                        end do
                    end do
                end if
            end do
        end
    ''', free_form = True)

    # apply
    k_loop = root_node.walk(Loop)[0]
    assert k_loop.variable.name == 'k'

    # patch
    with pytest.raises(Exception) as error:
        handle_kji_loop(k_loop, ['fx'])
    assert 'some k-use in the path' in str(error.value)

##########################################################
def test_handle_jki_loop():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            ! we will swap loop k<->i
            do j = 1, 100
                do k = 1, 10
                    do i = 1, 20
                        fx(i, j) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    j_loop = root_node.walk(Loop)[0]
    assert j_loop.variable.name == 'j'

    # patch
    ACCParallelTrans().apply(j_loop, options={})
    handle_jki_loop(j_loop, ['fx'], ['fx1d'])

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
  real, dimension(n,20) :: fx

  !$acc parallel default(present)
  do j = 1, 100, 1
    !$acc loop independent private(fx1d)
    do i = 1, 20, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

##########################################################
def test_handle_jik_loop_no_effect_single_loop():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            ! we will swap loop k<->i
            do j = 1, 10
                do i = 1, 10
                    do k = 1, 10
                        fx(i, j) = i+j+k
                    end do
                end do
            end do
        end
    ''', free_form = True)

    # apply
    j_loop = root_node.walk(Loop)[0]
    assert j_loop.variable.name == 'j'

    # patch
    ACCParallelTrans().apply(j_loop, options={})
    handle_jik_loop(j_loop, ['fx'], ['fx1d'])

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
  real, dimension(n,20) :: fx

  !$acc parallel default(present)
  do j = 1, 10, 1
    !$acc loop independent private(fx1d)
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

##########################################################
def test_handle_jik_many_child_loops():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            ! we will swap loop k<->i
            do j = 1, 10
              do i = 1, 10
                do k = 1, 10
                  fx(i, j) = i+j+k
                end do
              end do
            end do
            do j = 1, 10
              do i = 1, 10
                do k = 1, 10
                  fx(i, j) = i+j+k
                end do
              end do
            end do
        end
    ''', free_form = True)

    # apply
    j_loop = root_node.walk(Loop)[0]
    assert j_loop.variable.name == 'j'

    # patch
    ACCParallelTrans().apply(j_loop, options={})
    handle_jik_loop(j_loop, ['fx'], ['fx1d'])

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
  real, dimension(n,20) :: fx

  !$acc parallel default(present)
  do j = 1, 10, 1
    !$acc loop independent private(fx1d)
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

##########################################################
def test_handle_jik_many_child_loops_2():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source('''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fx(N, 20)
            implicit none
            ! we will swap loop k<->i
            do j = 1, 10
              do i = 1, 10
                do k = 1, 10
                  fx(i, j) = i+j+k
                end do
                do k = 1, 10
                  fx(i, j) = i+j+k
                end do
              end do
            end do
        end
    ''', free_form = True)

    # apply
    j_loop = root_node.walk(Loop)[0]
    assert j_loop.variable.name == 'j'

    # patch
    ACCParallelTrans().apply(j_loop, options={})
    handle_jik_loop(j_loop, ['fx'], ['fx1d'])

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
  real, dimension(n,20) :: fx

  !$acc parallel default(present)
  do j = 1, 10, 1
    !$acc loop independent private(fx1d)
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''


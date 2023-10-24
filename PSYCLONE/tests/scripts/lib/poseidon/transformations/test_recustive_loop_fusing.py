##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# from pytest package
import pytest
# psyclone
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.backend.fortran import FortranWriter
# internal
from scripts.lib.poseidon.transformations.recursive_loop_fusing import RecusiveLoopFusing, Loop, TransformationError

##########################################################
CODE_1='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
  enddo
  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
EXPECT_CODE_1='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_2='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
  enddo
  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo
  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff3
    enddo
  enddo

end subroutine test
'''

##########################################################
EXPECT_CODE_2='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
      vrhs(i,j) = cff2
      vrhs(i,j) = cff3
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_3='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
  enddo
  do j = 0, 20, 1
    do i = 5, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
EXPECTED_CODE_3='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
    do i = 5, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_4='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 10, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
  enddo
  do j = 0, 20, 1
    do i = 5, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_5='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1
    enddo
  enddo
  cff3 = 5
  do j = 0, 20, 1
    do i = 5, 10, 1
      vrhs(i,j) = cff2
    enddo
  enddo

end subroutine test
'''

##########################################################
def parse_apply_regen(code, trans):
    # load
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(code, free_form=True)

    # make trans
    top_loop = psyir_tree.walk(Loop, stop_type=Loop)
    trans.apply(top_loop)

    # write again
    writer = FortranWriter()
    final = writer(psyir_tree)

    # ok
    return final

##########################################################
def test_loop_fuse_basic():
    '''See if merge simple case'''
    # do
    assert parse_apply_regen(CODE_1, RecusiveLoopFusing()) == EXPECT_CODE_1

##########################################################
def test_loop_fuse_basic_2():
    '''See if merge simple case'''
    # do
    assert parse_apply_regen(CODE_2, RecusiveLoopFusing()) == EXPECT_CODE_2

##########################################################
def test_loop_fuse_partial_failure_1():
    '''See if merge simple case'''
    # do
    assert parse_apply_regen(CODE_3, RecusiveLoopFusing()) == EXPECTED_CODE_3

##########################################################
def test_loop_fuse_failure_3():
    '''See if merge simple case'''
    # do
    with pytest.raises(TransformationError):
        assert parse_apply_regen(CODE_4, RecusiveLoopFusing()) == CODE_4

##########################################################
def test_loop_fuse_failure_4():
    '''See if merge simple case'''
    # do
    with pytest.raises(TransformationError):
        assert parse_apply_regen(CODE_5, RecusiveLoopFusing()) == CODE_5

##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# from pytest package
import pytest
# psyclone
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.backend.fortran import FortranWriter
# internal
from scripts.poseidon.transformations.make_loop_single_top_level_operand import MakeLoopSingleTopLevelOperandTrans, Loop, TransformationError

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
      vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
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
      vrhs(i,j) = cff1 * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) + cff2 * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) + cff3 * vbar(i,j)
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
      vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) - cff3 * vbar(i,j)
    enddo
  enddo

end subroutine test
'''

##########################################################
EXPECT_CODE_3='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real :: cff1
  real :: cff2
  real :: cff3
  integer*4 :: i
  integer*4 :: j

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1 * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) + cff2 * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) - cff3 * vbar(i,j)
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

  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1 * vbar(i,j) * cff2 * vbar(i,j) * cff3 * vbar(i,j)
    enddo
  enddo

end subroutine test
'''

##########################################################
EXPECT_CODE_4='''subroutine test(cff1, cff2, cff3, vrhs, vbar)
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
      vrhs(i,j) = vrhs(i,j) * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) * cff2
      vrhs(i,j) = vrhs(i,j) * vbar(i,j)
      vrhs(i,j) = vrhs(i,j) * cff3
      vrhs(i,j) = vrhs(i,j) * vbar(i,j)
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
      cff1 = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
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
    top_loop = psyir_tree.walk(Loop)[0]
    trans.apply(top_loop)

    # write again
    writer = FortranWriter()
    final = writer(psyir_tree)

    # ok
    return final

##########################################################
def test_loop_single_op():
    '''If we apply on a loop with single line, should be not changed'''
    # do
    with pytest.raises(TransformationError) as error:
        assert parse_apply_regen(CODE_1, MakeLoopSingleTopLevelOperandTrans()) == CODE_1
    #assert "not contain binary operation" in str(error)

##########################################################
def test_loop_split_add():
    '''Shoud effectively split add ops'''
    # do
    assert parse_apply_regen(CODE_2, MakeLoopSingleTopLevelOperandTrans()) == EXPECT_CODE_2

##########################################################
def test_loop_split_add_sub():
    '''Shoud effectively split add/sub ops'''
    # do
    assert parse_apply_regen(CODE_3, MakeLoopSingleTopLevelOperandTrans()) == EXPECT_CODE_3

##########################################################
def test_loop_split_full_mul():
    '''Shoud effectively split add/sub ops'''
    # do
    assert parse_apply_regen(CODE_4, MakeLoopSingleTopLevelOperandTrans()) == EXPECT_CODE_4

##########################################################
def test_loop_not_split_scalar():
    '''If we apply on a loop with single line, should be not changed'''
    # do
    with pytest.raises(TransformationError) as error:
        assert parse_apply_regen(CODE_5, MakeLoopSingleTopLevelOperandTrans()) == CODE_5
    #assert "not contain binary operation" in str(error)

##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# from pytest package
import pytest
# psycone
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.backend.fortran import FortranWriter
# internal
from scripts.lib.poseidon.transformations.make_loop_single_assignment import MakeLoopSingleAssignmentTrans, Loop, TransformationError

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
      vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_2='''
subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
    REAL      :: vrhs(0:10, 0:20)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: dvom(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: om_v(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    do j = 0, 20, 1
        do i = 0, 10, 1
            vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
EXPECT_CODE_2='''subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real, dimension(0:10,0:20) :: dvom
  real, dimension(0:10,0:20) :: drhs
  real, dimension(0:10,0:20) :: om_v
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
  do j = 0, 20, 1
    do i = 0, 10, 1
      dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
    enddo
  enddo

end subroutine test
'''

##########################################################
CODE_3='''
subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
    REAL      :: vrhs(0:10, 0:20)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: dvom(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: om_v(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    do j = 0, 20, 1
        do i = 0, 10, 1
            cff1 = cff1 + cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
CODE_4='''
subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
    REAL      :: vrhs(0:10, 0:20)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: dvom(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: om_v(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    do j = 0, 20, 1
        do i = 0, 10, 1
            cff1 = cff1 + cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
        enddo
        do i = 0, 10, 1
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
CODE_5='''
subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
    REAL      :: vrhs(0:10, 0:20)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: dvom(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: om_v(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    cff1 = 5.0
    do j = 0, 20, 1
        if (cff1 == 3.0) then
            do i = 0, 10, 1
                cff1 = cff1 + cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
                dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
            enddo
        end if
    enddo
end
'''

##########################################################
CODE_6='''
subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
    REAL      :: vrhs(0:10, 0:20)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: dvom(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: om_v(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    do j = 0, 20, 1
        do i = 0, 10, 1
            vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
    do j = 0, 20, 1
        do i = 0, 10, 1
            vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
EXPECT_CODE_6='''subroutine test(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
  real, dimension(0:10,0:20) :: vrhs
  real, dimension(0:10,0:20) :: vbar
  real, dimension(0:10,0:20) :: dvom
  real, dimension(0:10,0:20) :: drhs
  real, dimension(0:10,0:20) :: om_v
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
  do j = 0, 20, 1
    do i = 0, 10, 1
      dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
    enddo
  enddo
  do j = 0, 20, 1
    do i = 0, 10, 1
      vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
      dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
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
def test_loop_single():
    '''If we apply on a loop with single line, should be not changed'''
    # do
    assert parse_apply_regen(CODE_1, MakeLoopSingleAssignmentTrans()) == CODE_1

##########################################################
def test_loop_two():
    '''Shoud effectively split the loops'''
    # do
    assert parse_apply_regen(CODE_2, MakeLoopSingleAssignmentTrans()) == EXPECT_CODE_2

##########################################################
def test_loop_scalar_in_loops():
    '''Detect error if scala assignement in loops (not supported yet)'''
    with pytest.raises(TransformationError) as error:
        parse_apply_regen(CODE_3, MakeLoopSingleAssignmentTrans())
    #assert "scalar assignement" in str(error)

##########################################################
def test_loop_complex_nesting():
    '''Check that we do not handle complex loop nesting'''
    with pytest.raises(TransformationError) as error:
        parse_apply_regen(CODE_4, MakeLoopSingleAssignmentTrans())

##########################################################
def test_loop_complex_nesting_with_if():
    '''Check that we do not handle complex if statement in loops'''
    with pytest.raises(TransformationError) as error:
        parse_apply_regen(CODE_5, MakeLoopSingleAssignmentTrans())

##########################################################
def test_loop_appy_in_middle_of_others():
    '''check that we insert at right place, if parent loop is followed by other ops'''
    assert parse_apply_regen(CODE_6, MakeLoopSingleAssignmentTrans()) == EXPECT_CODE_6

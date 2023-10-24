from poseidon_openacc_kernels_ifup import CrocoACCRaiseKernelThroughIf, trans
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from poseidon.dsl.helper import *
from psyclone.psyir.nodes.routine import Routine

def test_case_1():
    code = '''
    subroutine step2D_FB_tile(istr, jstr, jendr, zeta, knew)
        integer*4 :: istr
        integer*4 :: jstr
        integer*4 :: jendr
        integer*4 :: j
        integer*4 :: knew
        REAL :: zeta(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E, 4)
        REAL :: h(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E)
        REAL :: dnew(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E, 4)
        if (istr == 1) then
            do j = jstr - 1, jendr, 1
               dnew(istr - 1,j) = h(istr - 1,j) + zeta(istr - 1,j,knew)
            enddo
        end if
        return
    end subroutine step2D_FB_tile
    '''
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(code, free_form=False)

    kernels = extract_kernels_from_psyir(psyir_tree)
    kernels.make_acc_tranformation(False)
    kernels.merge_joinable_kernels()

    routines = psyir_tree.walk(Routine)
    for routine in routines:
        CrocoACCRaiseKernelThroughIf().apply(routine)
    
    assert FortranWriter()(psyir_tree) == '''subroutine step2D_FB_tile(istr, jstr, jendr, zeta, knew)
  integer*4 :: istr
  integer*4 :: jstr
  integer*4 :: jendr
  integer*4 :: knew
  REAL :: zeta(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E, 4)
  integer*4 :: j
  REAL :: h(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E)
  REAL :: dnew(0 : Lm + 1 + padd_X, 0 : Mm + 1 + padd_E, 4)

  !$acc kernels default(present) async(1)
  if (istr == 1) then
    do j = jstr - 1, jendr, 1
      dnew(istr - 1,j) = h(istr - 1,j) + zeta(istr - 1,j,knew)
    enddo
  end if
  !$acc end kernels
  return

end subroutine step2D_FB_tile
'''
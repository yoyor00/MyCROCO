##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
List the features which needs to be supported and needs some work to be supported.
'''

##########################################################
# from pytest package
import pytest
# psyclone
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.nodes import Routine, FileContainer, CodeBlock

##########################################################
CODE_COMMON = '''
      subroutine zetabc_tile(Istr,Iend,Jstr,Jend)
          implicit none
          integer*4 Istr,Iend,Jstr,Jend, i,j
          integer*4 IstrR,IendR,JstrR,JendR
          common /comm_setup_mpi1/ Lmmpi,Mmmpi
          if (istr.eq.1) then
              IstrR=Istr-1
          endif
      end
'''

##########################################################
CODE_PARAMETER = '''subroutine zetabc_tile()
  integer*4 :: size
  parameter (size=3)
  REAL :: h(0 : size + 1)


end subroutine zetabc_tile
'''

##########################################################
CODE_PARAMETER_REGEN = '''subroutine zetabc_tile()
  integer*4, parameter :: size = 3
  REAL :: h(0 : size + 1)


end subroutine zetabc_tile
'''

##########################################################
CODE_FIXED_FORM='''
      subroutine insert_node (lstr, node, nnodes, ierr)
      integer*4 lstr, ierr, i,j,k, lsffx, digits, power, ndots, idot(3)
     &                                        , node,  nnodes
      end
'''

##########################################################
CODE_WRITE_UNSUPPORTED='''
      subroutine step2d()
        implicit none
        write(stdout,'(A,F10.2)') ' VMAX (M/S) =   NaN'
      end
'''

##########################################################
CODE_STEP3D_ISSUE_1='''
      subroutine step3d_t (tile)

      integer tile
      real epsil
      parameter (epsil=1.E-16)

      call step3d_t_tile (epsil, tile)
      end
'''

##########################################################
def test_parse_fixed_form():
    '''
    When starting the free_form option was missing.
    '''
    fortran_reader = FortranReader()
    fortran_reader.psyir_from_source(CODE_FIXED_FORM, free_form = False)

##########################################################
@pytest.mark.xfail
def test_handle_common_keyword():
    '''
    Test if we at least accept the `common` keyword.
    '''
    fortran_reader = FortranReader()
    file_container = fortran_reader.psyir_from_source(CODE_COMMON, free_form = True)
    print(file_container.view())
    assert isinstance(file_container, FileContainer)
    subroutine = file_container.children[0]
    # Originally we have CodeBlock here because 'common' is not supported.
    assert not isinstance(block, CodeBlock)
    assert isinstance(subroutine, Routine)
    block = subroutine.children[0]

##########################################################
def test_handle_parameter_keyword():
    '''
    When starting the free_form option was missing.
    '''
    fortran_reader = FortranReader()
    file_container = fortran_reader.psyir_from_source(CODE_PARAMETER, free_form = True)
    writer = FortranWriter()
    result = writer(file_container)
    # we should have generated the same code.
    # Current issue : the parameter line is missing !
    assert result == CODE_PARAMETER_REGEN

##########################################################
def contains_block(node):
    if isinstance(node, CodeBlock):
        return True
    for child in node.children:
        if contains_block(child):
            return True
    return False

##########################################################
@pytest.mark.xfail
def test_handle_write_keyword():
    '''
    Test if write() call are still translated as CodeBlock
    '''
    fortran_reader = FortranReader()
    file_container = fortran_reader.psyir_from_source(CODE_WRITE_UNSUPPORTED, free_form = False)
    print(file_container.view())
    assert isinstance(file_container, FileContainer)
    subroutine = file_container.children[0]
    assert isinstance(subroutine, Routine)
    block = subroutine.children[0]
    # Originally we have CodeBlock here because 'common' is not supported.
    assert not contains_block(block)

##########################################################
def test_parse_step_3d_1():
    '''
    When starting the free_form option was missing.
    '''
    fortran_reader = FortranReader()
    fortran_reader.psyir_from_source(CODE_STEP3D_ISSUE_1, free_form = False)

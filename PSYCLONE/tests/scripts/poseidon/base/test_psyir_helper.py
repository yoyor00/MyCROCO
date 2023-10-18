##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# psyclone
from psyclone.psyir.frontend.fortran import FortranReader

CODE='''
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

def test_extract_var_ref_list():
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(CODE)
    res = psyir_helpers.extract_var_ref_list(psyir_tree, 0)
    assert res == [
         {'mode': 'W', 'var': 'vrhs'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'cff1'},
         {'mode': 'R', 'var': 'vbar'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'cff2'},
         {'mode': 'R', 'var': 'vbar'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'cff3'},
         {'mode': 'R', 'var': 'vbar'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'W', 'var': 'dvom'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'drhs'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'drhs'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'om_v'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
         {'mode': 'R', 'var': 'vrhs'},
         {'mode': 'R', 'var': 'i'},
         {'mode': 'R', 'var': 'j'},
    ]

def test_extract_var_ref_list_reduced():
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(CODE)
    vars = psyir_helpers.extract_var_ref_list_reduced(psyir_tree)
    assert vars == {
        'cff1': 'R',
        'cff2': 'R',
        'cff3': 'R',
        'drhs': 'R',
        'dvom': 'W',
        'i': 'R',
        'j': 'R',
        'om_v': 'R',
        'vbar': 'R',
        'vrhs': 'RW'
    }

def test_is_using_var():
    reader = FortranReader()
    psyir_tree = reader.psyir_from_source(CODE)
    assert psyir_helpers.is_using_var(psyir_tree, 'cff1') == True
    assert psyir_helpers.is_using_var(psyir_tree, 'i') == True
    assert psyir_helpers.is_using_var(psyir_tree, 'k') == False

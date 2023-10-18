##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper
functions.
'''

##########################################################
# pytest package
import pytest
# psyclone
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
# internal
from scripts.extensions.loops import *

##########################################################
VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

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
    handle_jki_loop(j_loop, ['fx'])

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

  do j = 1, 100, 1
    do i = 1, 20, 1
      do k = 1, 10, 1
        fx(i,j) = i + j + k
      enddo
    enddo
  enddo

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
    handle_jik_loop(j_loop, ['fx'])

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

  do j = 1, 10, 1
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo

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
    handle_jik_loop(j_loop, ['fx'])

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

  do j = 1, 10, 1
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo

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
    handle_jik_loop(j_loop, ['fx'])

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

  do j = 1, 10, 1
    do i = 1, 10, 1
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
      do k = 1, 10, 1
        fx1d(j) = i + j + k
      enddo
    enddo
  enddo

end subroutine step3d_t
'''

##########################################################
def helper_load_snippet(kind: str, step: str, name: str, generated_node: Node = None, dump: bool = False) -> str:
    # calc path
    test_dir = os.path.dirname(os.path.realpath(__file__))
    snipped_dir = os.path.join(test_dir, "croco-loops-snippets/", step, kind)
    fname = os.path.join(snipped_dir, f"snippet-{name}.F")

    # check it matches
    if generated_node != None:
        generated_source = FortranWriter()(generated_node)
    
        # dump if enabled
        if dump:
            os.makedirs(snipped_dir, exist_ok=True)
            with open(fname, 'w+') as fp:
                fp.write(generated_source)

    # load
    with open(fname, 'r') as fp:
        current_source = fp.read()

    # check it matches
    if generated_node != None:
        assert current_source == generated_source

    # ok
    return current_source

##########################################################
def helper_gen_var_decl(vars: dict, common: bool = False) -> str:
    # to aggregate
    decl = []

    # loop on scalars
    scalars = vars['scalars_int']
    for vname in scalars:
        decl.append(f'integer*4 {vname}')

    # loop on scalars
    scalars = vars['scalars']
    for vname in scalars:
        decl.append(f'real*4 {vname}')

    # loop on 1d arrays
    scalars = vars['1d']
    for vname in scalars:
        decl.append(f'real*4 {vname}(0:10)')

    # loop on 2d arrays
    scalars = vars['2d']
    for vname in scalars:
        decl.append(f'real*4 {vname}(0:10,0:10)')

    # loop on 3d arrays
    scalars = vars['3d']
    for vname in scalars:
        decl.append(f'real*4 {vname}(0:10,0:10,0:10)')

    # loop on 4d arrays
    scalars = vars['4d']
    for vname in scalars:
        decl.append(f'real*4 {vname}(0:10,0:10,0:10,0:10)')

    # loop on 5d arrays
    scalars = vars['5d']
    for vname in scalars:
        decl.append(f'real*4 {vname}(0:10,0:10,0:10,0:10,0:10)')
  
    # common
    if common:
        for dim in vars:
            for vname in vars[dim]:
                decl.append(f'common {vname}')

    # ok
    return '\n'.join(decl)

##########################################################
LOOP_USED_VARS={
    'scalars_int': [
        'itrc', 'nt', 'i', 'j', 'k', 'nadv', 'iend', 'lm', 'n', 'jstr', 'jend',
        'istr', 'mm', 'iic', 'ntstart', 'nnew', 'nstp',
        'indx', 'jstrv', 'istru', 'itemp',
        'nrhs',
        'jendr', 'istrr', 'iendr', 'knew', 'jstrr'],
    'scalars': [
        'cff', 'dt', 'gamma',
        'cff1', 'cff2', 'cdt', 'epsil',
        'cfr', 'g', 'grho', 'onefifth', 'halfgrho', 'onetwelfth',
        'aa', 'cc'],
    '1d': ['fc1d', 'cf1d', 'dc1d', 'bc1d', 'dz1d', 'dr1d'],
    '2d': [
        'fx', 'work', 'fe', 'pm', 'pn', 'fc', 'cf', 'dz', 'dr', 'on_u',
        'rufrc', 'rvfrc', 'bc', 'dc', 'om_v', 'dv_avg2', 'du_avg2'
    ],
    '3d': [
        'huon', 'hvom', 'hz', 'hz_half', 'rv', 'we', 'stflx',
        'btflx', 'z_r', 'vfx_3d', 'vfe_3d', 'ufx_3d', 'ufe_3d', 'rho',
        'p', 'z_w', 'du_avg1', 'ubar', 'dv_avg1',
        'vbar', 'romega', 'ru', 'work_3d', 'hz_bak', 'fx_3d', 'fe_3d',
        'akv'
    ],
    '4d': ['v', 'u', 'akt'],
    '5d': ['t']
}
LOOP_PARAMETERS = [
    # kji loops
    ('kji', 'pre_step3d_tile-1129'),
    ('kji', 'pre_step3d_tile-1236'),
    ('kji', 'pre_step3d_tile-1491'),
    ('kji', 'pre_step3d_tile-115'),
    ('kji', 'pre_step3d_tile-1384'),
    ('kji', 'pre_step3d_tile-264'),

    # jik loops
    ('jik', 'omega_tile-158'),
    ('jik', 'step3d_uv2_tile-399'),
    ('jik', 'step3d_uv2_tile-73'),
    ('jik', 'rhs3d_tile-2289'),
    ('jik', 'step3d_uv2_tile-737'),
    ('jik', 'step3d_uv2_tile-926'),

    # jki loops
    ('jki', 'pre_step3d_tile-873'),
    ('jki', 'prsgrd_tile-83'),
    ('jki', 'rhs3d_tile-1781'),
    ('jki', 'rhs3d_tile-2086'),
    ('jki', 'step3d_t_tile-526')
]
@pytest.mark.parametrize(("type", "snippet_name"), LOOP_PARAMETERS)
def test_handle_loops_from_croco(type: str, snippet_name: str):
    '''
    Check whether it inserts the vars in subroutine.
    '''

    # load snippet
    snippet = helper_load_snippet(f'{type}-loops', 'origin', snippet_name)

    # embed in routine
    full_source = f'''\
subroutine snippet()
implicit none
{helper_gen_var_decl(LOOP_USED_VARS)}
{snippet}
end
'''

    # debug
    print(full_source)

    # parse to get IR tree
    root_node: Node
    root_node = FortranReader().psyir_from_source(full_source, free_form = True)

    # To use when first importing a new snippet so it generte the files
    # CAUTION : need to check that the expected sources are correct as it is the reference !
    # DO NOT FORGET : to make it again false after generating the expectation otherwise it makes the test useless !
    dump = False

    # ensuite we are sync with the prepared sources so we can easily kompare on the dirs
    prep_source = helper_load_snippet(f'{type}-loops', 'prepared', snippet_name, generated_node = root_node, dump = dump)

    # patch
    if type == 'kji':
        k_loop = get_first_loop_on(root_node.walk(Loop)[0], 'k')
        assert k_loop.variable.name == 'k'
        handle_kji_loop(k_loop, VARS_3D)
    elif type == 'jki':
        j_loop = get_first_loop_on(root_node.walk(Loop)[0], 'j')
        assert j_loop.variable.name == 'j'
        handle_jki_loop(j_loop, VARS_1D)
    elif type == 'jik':
        j_loop = get_first_loop_on(root_node.walk(Loop)[0], 'j')
        assert j_loop.variable.name == 'j'
        handle_jik_loop(j_loop, VARS_1D, do_k_loop_fuse=True)
    else:
        assert False

    # regen
    gen_source = FortranWriter()(root_node)
    #print(gen_source)

    # load expectation
    expected_source = helper_load_snippet(f'{type}-loops', 'expected', snippet_name, generated_node = root_node, dump = dump)

    # check
    assert gen_source == expected_source

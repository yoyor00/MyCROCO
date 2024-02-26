##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper
functions.
'''

##########################################################
# python
import os
import pytest
import subprocess
import tempfile
import shutil
# Psyclone
from psyclone.psyir.nodes import Node, Loop, Routine, Call
# internal
from scripts.lib.extensions import acc
from scripts.lib.extensions.loops import get_first_loop_on, handle_kji_loop, handle_jki_loop, handle_jik_loop, dump_named_node_as_source
from scripts.lib.croco import psyacc, step3d
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
# internal poseidon
from scripts.lib.poseidon.dsl.helper import extract_kernels_from_psyir
# local test dir
from .test_loops import helper_load_snippet, helper_gen_var_decl

##########################################################
VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

##########################################################
def helper_gen_init(vars: dict) -> str:
    # to aggregate
    decl = []

    decl.append('integer, save :: seed_size')
    decl.append('integer, save, dimension(:), allocatable :: seeder')
    decl.append('logical, save :: first_call = .TRUE.')

    # random
    #decl.append(f'call random_init(.true., .true.)')
    decl.append(f'if (first_call) then')
    decl.append(f'  call random_seed(size=seed_size)')
    decl.append(f'  allocate(seeder(seed_size))')
    decl.append(f'  call random_seed(get=seeder)')
    decl.append(f'  first_call = .FALSE.')
    decl.append(f'else')
    decl.append(f'  call random_seed(put=seeder)')
    decl.append(f'endif')

    # loop on scalars
    scalars = vars['scalars']
    #for vname in scalars:
    #    decl.append(f'integer*4 {vname}')
    decl.append(f'n = 10')
    decl.append(f'jstr = 2')
    decl.append(f'jend = 9')
    decl.append(f'istr = 2')
    decl.append(f'iend = 9')
    decl.append(f'jmin = jstr-2')
    decl.append(f'jmax = jend+1')
    decl.append(f'imin = istr-2')
    decl.append(f'imax = iend+1')

    # loop on 1d arrays
    scalars = vars['1d']
    for vname in scalars:
        decl.append(f'call random_number({vname})')

    # loop on 2d arrays
    scalars = vars['2d']
    for vname in scalars:
        decl.append(f'call random_number({vname})')

    # loop on 3d arrays
    scalars = vars['3d']
    for vname in scalars:
        decl.append(f'call random_number({vname})')

    # loop on 4d arrays
    scalars = vars['4d']
    for vname in scalars:
        decl.append(f'call random_number({vname})')

    # loop on 5d arrays
    scalars = vars['5d']
    for vname in scalars:
        decl.append(f'call random_number({vname})')

    # ok
    return '\n'.join(decl)

##########################################################
def helper_gen_var_decl_saved(vars: dict) -> str:
    # rename cf => saved_cf
    renamed = {}
    for cat, names in vars.items():
        new_names = []
        for name in names:
            new_names.append(f"saved_{name}")
        renamed[cat] = new_names

    # generate vars
    return helper_gen_var_decl(renamed, common=True)

##########################################################
def helper_gen_var_saving(vars: dict) -> str:
    # to aggregate
    decl = []

    # loop on all vars
    for cat, names in vars.items():
        for name in names:
            decl.append(f"saved_{name} = {name}")

     # ok
    return '\n'.join(decl)

##########################################################
def helper_gen_check(vars: dict, skip_check: dict) -> str:
    # to aggregate
    decl = []

    # loop on scalars
    scalars = vars['scalars']
    #for vname in scalars:
    #    decl.append(f'integer*4 {vname}')
    decl.append(f'n = 10')
    decl.append(f'jstr = 2')
    decl.append(f'jend = 9')
    decl.append(f'istr = 2')
    decl.append(f'iend = 9')

    # loop on 1d arrays
    scalars = vars['1d']
    #for vname in scalars:
    #    decl.append(f'  random_number({vname})')

    # loop on 2d arrays
    for dim in range(1,5):
        scalars = vars[f'{dim}d']
        for vname in scalars:
            if skip_check and not vname in skip_check:
                if not vname.endswith('_3d') and not vname.endswith('1d'):
                    decl.append(f'if (all({vname} .ne. saved_{vname})) then')
                    decl.append(f"  write(*,*) 'Invalid value in {vname} !'")
                    decl.append(f'  call exit(1)')
                    decl.append(f'endif')

    # ok
    return '\n'.join(decl)

##########################################################
LOOP_USED_VARS={
    'scalars_int': [
        'itrc', 'nt', 'i', 'j', 'k', 'nadv', 'iend', 'lm', 'n', 'jstr', 'jend',
        'istr', 'mm', 'iic', 'ntstart', 'nnew', 'nstp',
        'indx', 'jstrv', 'istru', 'itemp',
        'nrhs',
        'jendr', 'istrr', 'iendr', 'knew', 'jstrr',
        'jmin', 'jmax', 'imin', 'imax',
        'kstp', 'kbak', 'kold', 'levsfrc',
        'levbfrc'],
    'scalars': [
        'cff', 'dt', 'gamma',
        'cff1', 'cff2', 'cff3', 'cdt', 'epsil',
        'cfr', 'g', 'grho', 'onefifth', 'halfgrho', 'onetwelfth',
        'aa', 'cc',
        'rdrg', 'rdrg2',
        'dunew', 'dvnew',
        'vmax', 'dtfast', 'cff0', 'gamma2', 'cffx', 'curvx', 'cffe', 'curve'],
    '1d': ['fc1d', 'cf1d', 'dc1d', 'bc1d', 'dz1d', 'dr1d'],
    '2d': [
        'fx', 'work', 'fe', 'pm', 'pn', 'fc', 'cf', 'dz', 'dr', 'on_u',
        'rufrc', 'rvfrc', 'bc', 'dc', 'om_v', 'dv_avg2', 'du_avg2',
        'rubar', 'ufx', 'ufe', 'rvbar', 'vfe', 'vfx', 'drhs', 'vrhs', 'fomn',
        'urhs', 'om_u', 'on_v', 'rhos', 'rhoa', 'zeta_new', 'h', 'duon', 'dnew',
        'dvom', 'zt_avg1', 'om_r', 'on_r', 'pm_u', 'pn_u', 'pm_v', 'pn_v',
        'sustr', 'svstr', 'bustr', 'bvstr'
    ],
    '3d': [
        'huon', 'hvom', 'hz', 'hz_half', 'rv', 'we', 'stflx',
        'btflx', 'z_r', 'vfx_3d', 'vfe_3d', 'ufx_3d', 'ufe_3d', 'rho',
        'p', 'z_w', 'du_avg1', 'ubar', 'dv_avg1',
        'vbar', 'romega', 'ru', 'work_3d', 'hz_bak', 'fx_3d', 'fe_3d',
        'akv', 'rufrc_bak', 'rvfrc_bak', 'zeta', 'wrk1_3d', 'wrk2_3d'
    ],
    '4d': ['v', 'u', 'akt'],
    '5d': ['t']
}
LOOP_PARAMETERS = [
    # kji loops
    ('kji', 'pre_step3d_tile-1129', []),
    ('kji', 'pre_step3d_tile-1236', []),
    ('kji', 'pre_step3d_tile-1491', []),
    ('kji', 'pre_step3d_tile-115', []),
    ('kji', 'pre_step3d_tile-1384', []),
    ('kji', 'pre_step3d_tile-264', []),

    # jik loops
    ('jik', 'omega_tile-158', []),
    ('jik', 'step3d_uv2_tile-399', []),
    ('jik', 'step3d_uv2_tile-73', []),
    ('jik', 'rhs3d_tile-2289', []),
    ('jik', 'step3d_uv2_tile-737', []),
    ('jik', 'step3d_uv2_tile-926', []),

    # jki loops
    ('jki', 'pre_step3d_tile-873', []),
    ('jki', 'prsgrd_tile-83', ['dz']),
    ('jki', 'rhs3d_tile-1781', []),
    ('jki', 'rhs3d_tile-2086', []),
    ('jki', 'step3d_t_tile-526', []),

    # step2d kernel openacc loops
    ('kernel-openacc', 'step2D_FB_tile-1054', []),
    ('kernel-openacc', 'step2D_FB_tile-1122', []),
    ('kernel-openacc', 'step2D_FB_tile-1179', []),
    ('kernel-openacc', 'step2D_FB_tile-1212', []),
    ('kernel-openacc', 'step2D_FB_tile-1252', []),
    ('kernel-openacc', 'step2D_FB_tile-1336', []),
    ('kernel-openacc', 'step2D_FB_tile-1428', []),
    ('kernel-openacc', 'step2D_FB_tile-1464', []),
    ('kernel-openacc', 'step2D_FB_tile-1556', []),
    ('kernel-openacc', 'step2D_FB_tile-1608', []),
    ('kernel-openacc', 'step2D_FB_tile-1660', []),
    ('kernel-openacc', 'step2D_FB_tile-1741', []),
    ('kernel-openacc', 'step2D_FB_tile-185', []),
    ('kernel-openacc', 'step2D_FB_tile-1927', []),
    ('kernel-openacc', 'step2D_FB_tile-1971', []),
    ('kernel-openacc', 'step2D_FB_tile-2064', []),
    ('kernel-openacc', 'step2D_FB_tile-2176', []),
    ('kernel-openacc', 'step2D_FB_tile-2212', []),
    ('kernel-openacc', 'step2D_FB_tile-2248', []),
    ('kernel-openacc', 'step2D_FB_tile-227', []),
    ('kernel-openacc', 'step2D_FB_tile-2284', []),
    ('kernel-openacc', 'step2D_FB_tile-2327', []),
    ('kernel-openacc', 'step2D_FB_tile-2376', []),
    ('kernel-openacc', 'step2D_FB_tile-2434', []),
    ('kernel-openacc', 'step2D_FB_tile-2481', []),
    ('kernel-openacc', 'step2D_FB_tile-2539', []),
    ('kernel-openacc', 'step2D_FB_tile-2590', []),
    ('kernel-openacc', 'step2D_FB_tile-2646', []),
    ('kernel-openacc', 'step2D_FB_tile-2697', []),
    ('kernel-openacc', 'step2D_FB_tile-2752', []),
    ('kernel-openacc', 'step2D_FB_tile-2773', []),
    ('kernel-openacc', 'step2D_FB_tile-294', []),
    ('kernel-openacc', 'step2D_FB_tile-414', []),
    ('kernel-openacc', 'step2D_FB_tile-470', []),
    ('kernel-openacc', 'step2D_FB_tile-557', []),
    ('kernel-openacc', 'step2D_FB_tile-616', []),
    ('kernel-openacc', 'step2D_FB_tile-673', []),
    ('kernel-openacc', 'step2D_FB_tile-734', []),
    ('kernel-openacc', 'step2D_FB_tile-910', []),
    ('kernel-openacc', 'step2D_FB_tile-982', []),
    ('kernel-openacc', 'u2dbc_tile-102', []),
    ('kernel-openacc', 'u2dbc_tile-126', []),
    ('kernel-openacc', 'u2dbc_tile-155', []),
    ('kernel-openacc', 'u2dbc_tile-80', []),
    ('kernel-openacc', 'v2dbc_tile-102', []),
    ('kernel-openacc', 'v2dbc_tile-126', []),
    ('kernel-openacc', 'v2dbc_tile-155', []),
    ('kernel-openacc', 'v2dbc_tile-80', []),
    ('kernel-openacc', 'zetabc_tile-107', []),
    ('kernel-openacc', 'zetabc_tile-136', []),
    ('kernel-openacc', 'zetabc_tile-163', []),
    ('kernel-openacc', 'zetabc_tile-80', []),

    # step3d specific
    ('kernel-step3d', 'omega_tile-158', []),
    ('kernel-step3d', 'omega_tile-294', []),
    ('kernel-step3d', 'omega_tile-320', []),
    ('kernel-step3d', 'omega_tile-348', []),
    ('kernel-step3d', 'omega_tile-374', []),
    ('kernel-step3d', 'omega_tile-406', []),
    ('kernel-step3d', 'omega_tile-431', []),
    ('kernel-step3d', 'omega_tile-458', []),
    ('kernel-step3d', 'omega_tile-485', []),
    ('kernel-step3d', 'omega_tile-73', []),
    ('kernel-step3d', 'omega_tile-98', []),
    ('kernel-step3d', 'pre_step3d_tile-1072', []),
    ('kernel-step3d', 'pre_step3d_tile-1100', []),
    ('kernel-step3d', 'pre_step3d_tile-115', []),
    ('kernel-step3d', 'pre_step3d_tile-1206', []),
    ('kernel-step3d', 'pre_step3d_tile-1325', []),
    ('kernel-step3d', 'pre_step3d_tile-1353', []),
    ('kernel-step3d', 'pre_step3d_tile-1459', []),
    #('kernel-step3d', 'pre_step3d_tile-1578', []),
    ('kernel-step3d', 'pre_step3d_tile-1603', []),
    ('kernel-step3d', 'pre_step3d_tile-264', []),
    ('kernel-step3d', 'pre_step3d_tile-788', []),
    ('kernel-step3d', 'prsgrd_tile-449', []),
    ('kernel-step3d', 'prsgrd_tile-489', []),
    ('kernel-step3d', 'prsgrd_tile-83', []),
    ('kernel-step3d', 'rhs3d_tile-128', []),
    ('kernel-step3d', 'rhs3d_tile-1679', []),
    ('kernel-step3d', 'rhs3d_tile-1730', []),
    ('kernel-step3d', 'rhs3d_tile-1781', []),
    ('kernel-step3d', 'rhs3d_tile-179', []),
    ('kernel-step3d', 'rhs3d_tile-2041', []),
    ('kernel-step3d', 'rhs3d_tile-2092', []),
    ('kernel-step3d', 'rhs3d_tile-2143', []),
    ('kernel-step3d', 'rhs3d_tile-230', []),
    ('kernel-step3d', 'rhs3d_tile-2403', []),
    ('kernel-step3d', 'rhs3d_tile-309', []),
    ('kernel-step3d', 'rhs3d_tile-330', []),
    ('kernel-step3d', 'rhs3d_tile-363', []),
    ('kernel-step3d', 'rhs3d_tile-414', []),
    ('kernel-step3d', 'rhs3d_tile-465', []),
    ('kernel-step3d', 'rhs3d_tile-544', []),
    ('kernel-step3d', 'rhs3d_tile-73', []),
    ('kernel-step3d', 'rhs3d_tile-94', []),
    #('kernel-step3d', 'step3d_t_tile-1082', []),
    ('kernel-step3d', 'step3d_t_tile-1095', []),
    ('kernel-step3d', 'step3d_t_tile-1102', []),
    ('kernel-step3d', 'step3d_t_tile-500', []),
    ('kernel-step3d', 'step3d_t_tile-76', []),
    ('kernel-step3d', 'step3d_uv2_tile-1203', []),
    ('kernel-step3d', 'step3d_uv2_tile-506', []),
    ('kernel-step3d', 'step3d_uv2_tile-73', []),
    ('kernel-step3d', 'step3d_uv2_tile-951', [])
]
#LOOP_PARAMETERS = [
#    ('kji', 'pre_step3d_tile-1129'),
#]
@pytest.mark.parametrize(("type", "snippet_name", "skip_check"), LOOP_PARAMETERS)
def test_running_reshaped_kernels(type: str, snippet_name: str, skip_check: list):
    # load snippet
    snippet = helper_load_snippet(f'{type}-loops', 'origin', snippet_name)
    snippet_func = helper_load_snippet(f'functions', 'extra', 't3dbc_tile')

    # embed in routine
    full_source = f'''\
program test_snippet
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {helper_gen_var_decl_saved(LOOP_USED_VARS)}

    call init_vars()
    call snippet_cpu_default_code()
    call save_vars()

    call init_vars()
    call snippet_reshaped_code()
    call post_check()
end

subroutine init_vars()
    implicit none
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {helper_gen_init(LOOP_USED_VARS)}
end

subroutine post_check()
    implicit none
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {helper_gen_var_decl_saved(LOOP_USED_VARS)}
    {helper_gen_check(LOOP_USED_VARS, skip_check)}
end

subroutine save_vars()
    implicit none
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {helper_gen_var_decl_saved(LOOP_USED_VARS)}
    {helper_gen_var_saving(LOOP_USED_VARS)}
end

subroutine snippet_cpu_default_code()
    implicit none
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {snippet}
end

subroutine snippet_reshaped_code()
    implicit none
    {helper_gen_var_decl(LOOP_USED_VARS, common=True)}
    {snippet}
end
'''

    # debug
    #print(full_source)

    # dump
    with open(os.path.expanduser("~/source.orig.F90"), "w+") as fp:
        fp.write(full_source)

    # parse to get IR tree
    top_root_node: Node
    top_root_node = FortranReader().psyir_from_source(full_source, free_form = True)

    # get sub function
    root_node: Node
    for routine in top_root_node.walk(Routine):
        if routine.name == 'snippet_reshaped_code':
            root_node = routine
            break

    # dump before
    dump_named_node_as_source(snippet_name, routine, f"{type}-orig", in_dir=os.path.expanduser("~/snippets/"))

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
    elif type == 'kernel-openacc':
        assert isinstance(root_node, Routine)
        routine: Routine = root_node
        psyacc.apply_psyclone_original_trans(root_node)
    elif type == 'kernel-step3d':
        assert isinstance(root_node, Routine)
        routine: Routine = root_node
        step3d.apply_step3d_routine_trans(root_node, None)
    else:
        assert False

    # dump before
    dump_named_node_as_source(snippet_name, routine, f"{type}-transformed", in_dir=os.path.expanduser("~/snippets/"))

    # inject acc kernel directives
    #kernels = extract_kernels_from_psyir(top_root_node)
    #kernels.make_acc_tranformation(False)
    #kernels.merge_joinable_kernels()

    # regen
    gen_source = FortranWriter()(top_root_node)
    #print(gen_source)

    # save & compiler & run
    exename = '/tmp/try_compiler_snippet'
    with tempfile.NamedTemporaryFile("w+", suffix=".F90") as fp_source:
        # write
        fp_source.write(gen_source)
        fp_source.flush()

        # compiler
        try:
            subprocess.run(['gfortran', fp_source.name, '-o', exename, '-ffree-line-length-none', '-Wall', '-Werror', '-Wno-unused-variable'], check=True)

            # run
            subprocess.run([exename], check=True)
        except Exception    as e:
            #subprocess.run(["cat", fp_source.name], check=True)
            shutil.copyfile(fp_source.name, os.path.expanduser("~/source.F90"))
            raise e

        # compiler
        #subprocess.run(['nvfortran', fp.name, '-o', exename, '-Wall', '-Werror', '-acc=gpu'], check=True)

        # run
        #subprocess.run([exename], check=True)

##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some basic unit test check check the transformation device helper functions.
'''

##########################################################
from tempfile import NamedTemporaryFile
from ..kernels import *
from psyclone.configuration import Config
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.transformations import ACCParallelTrans
from psyclone.parse.algorithm import parse
from psyclone.psyGen import PSyFactory
from poseidon.dsl.helper import extract_kernels_from_psyir

##########################################################
def helper_reproduce_psyclone_psy(source: str):
    # dump in config file
    with NamedTemporaryFile("w+") as fp:
        # dump
        fp.write(source)
        fp.flush()

        # get psyclone config
        api = "nemo"
        distributed_memory = Config.get().distributed_memory
        ast, invoke_info = parse(fp.name, api=api, invoke_name="invoke", kernel_paths=False, line_length=80)
        psy = PSyFactory(api, distributed_memory=distributed_memory).create(invoke_info)

    # ok
    return psy

##########################################################
def helper_get_psy_and_kernels(source: str) -> (object, KernelList):
    psy = helper_reproduce_psyclone_psy(source)

    # apply parallel
    ACCParallelTrans().apply(psy.container.walk(Loop)[0].parent)

    # extranct kernels
    kernels = extract_kernels_from_psyir(psy.container)

    # ok
    return psy, kernels

##########################################################
def test_apply_acc_loop_and_collapse():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    source = '''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real tt(N, 10, 20)
            implicit none
            do k = 1, 10
                do j = 1, 10
                    do i = 1, 20
                        tt(i, j, k) = i+j+k
                    end do
                end do
            end do
            do j = 1, 10
                do i = 1, 20
                    tt(i, j, k) = i+j+k
                end do
            end do
        end
    '''

    # get kernels
    psy, kernels = helper_get_psy_and_kernels(source)

    # apply
    apply_acc_loop_and_collapse(kernels, options={})

    # regen
    gen_source = FortranWriter()(psy.container)
    print(gen_source)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  real, dimension(n,10,20) :: tt

  !$acc parallel default(present)
  !$acc loop independent collapse(3)
  do k = 1, 10, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        tt(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$acc loop independent collapse(2)
  do j = 1, 10, 1
    do i = 1, 20, 1
      tt(i,j,k) = i + j + k
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

##########################################################
def test_apply_acc_fetch_vars():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    # parse to get IR tree
    source = '''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real fc1d(N, 10, 20)
            implicit none
            do k = 1, 100
                do j = 1, 10
                    do i = 1, 20
                        fc1d(i, j, k) = i+j+k
                    end do
                end do
            end do
        end
    '''

    # get kernels
    psy, kernels = helper_get_psy_and_kernels(source)

    # patch
    apply_acc_fetch_vars(psy)

    # regen
    gen_source = FortranWriter()(psy.container)
    print(gen_source)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  real, dimension(n,10,20) :: fc1d

  !$acc enter data copyin(fc1d,i,j,k)
  !$acc parallel default(present)
  do k = 1, 100, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        fc1d(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$acc end parallel

end subroutine step3d_t
'''

##########################################################
def test_apply_acc_kernels():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    source = '''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real tt(N, 10, 20)
            implicit none
            do k = 1, 10
                do j = 1, 10
                    do i = 1, 20
                        tt(i, j, k) = i+j+k
                    end do
                end do
            end do
            do j = 1, 10
                do i = 1, 20
                    tt(i, j, k) = i+j+k
                end do
            end do
        end
    '''

    # get kernels
    psy = helper_reproduce_psyclone_psy(source)

    # apply
    collapse = False
    joinable = False
    options = {'independent': False}
    routine = psy.container.walk(Routine)[0]
    apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable, options=options)

    # regen
    gen_source = FortranWriter()(psy.container)
    print(gen_source)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  real, dimension(n,10,20) :: tt

  !$acc kernels default(present)
  do k = 1, 10, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        tt(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$acc end kernels
  !$acc kernels default(present)
  do j = 1, 10, 1
    do i = 1, 20, 1
      tt(i,j,k) = i + j + k
    enddo
  enddo
  !$acc end kernels

end subroutine step3d_t
'''

##########################################################
def test_apply_acc_kernels():
    '''
    Check whether it inserts the vars in subroutine.
    '''
    source = '''
        subroutine step3d_t(N)
            integer*4 N, i, j, k
            real tt(N, 10, 20)
            implicit none
            do k = 1, 10
                do j = 1, 10
                    do i = 1, 20
                        tt(i, j, k) = i+j+k
                    end do
                end do
            end do
            do j = 1, 10
                do i = 1, 20
                    tt(i, j, k) = i+j+k
                end do
            end do
        end
    '''

    # get kernels
    psy = helper_reproduce_psyclone_psy(source)

    # apply
    collapse = True
    joinable = True
    options = {'independent': False}
    routine = psy.container.walk(Routine)[0]
    apply_acc_kernel(routine, collapse, ignore_loops=['itrc'], merge_joinable=joinable, options=options)

    # regen
    gen_source = FortranWriter()(psy.container)
    print(gen_source)

    # check
    assert gen_source == '''\
subroutine step3d_t(n)
  integer*4 :: n
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  real, dimension(n,10,20) :: tt

  !$acc kernels default(present)
  !$acc loop collapse(3)
  do k = 1, 10, 1
    do j = 1, 10, 1
      do i = 1, 20, 1
        tt(i,j,k) = i + j + k
      enddo
    enddo
  enddo
  !$acc loop collapse(2)
  do j = 1, 10, 1
    do i = 1, 20, 1
      tt(i,j,k) = i + j + k
    enddo
  enddo
  !$acc end kernels

end subroutine step3d_t
'''

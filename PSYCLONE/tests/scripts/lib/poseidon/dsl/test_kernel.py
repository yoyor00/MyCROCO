##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# internal
from scripts.lib.poseidon.dsl.helper import *
from scripts.lib.poseidon.dsl.kernel_pre_scheduler import KernelPreScheduler

##########################################################
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

    ! kernel 1
    do j = 0, 20, 1
        do i = 0, 10, 1
            vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo

    ! kernel 2
    do j = 0, 20, 1
        do i = 0, 10, 1
            vrhs(i,j) = vrhs(i,j) + 2.0d0 * vbar(i,j)
            dvom(i,j) = dvom(i,j) + 0.5d0 * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
KERNEL_1='''================= KERNEL ==================
IN    : cff1, cff2, cff3, drhs, dvom, om_v, vbar, vrhs
OUT   : dvom, vrhs
MASKS : (i, j - 1), (i, j)
------------------ ORIG -------------------
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
    dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo

-------------------------------------------
'''

##########################################################
KERNEL_2='''================= KERNEL ==================
IN    : cff1, cff2, cff3, drhs, dvom, om_v, vbar, vrhs
OUT   : dvom, vrhs
MASKS : (i, j - 1), (i, j)
------------------ ORIG -------------------
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
    dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo

-------------------------------------------
'''

##########################################################
KERNEL_3='''================= KERNEL ==================
IN    : cff1, cff2, cff3, drhs, dvom, om_v, vbar, vrhs
OUT   : dvom, vrhs
MASKS : (i, j - 1), (i, j)
------------------ ORIG -------------------
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
    dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo

-------------------------------------------
================= KERNEL ==================
IN    : dvom, vbar, vrhs
OUT   : dvom, vrhs
MASKS : (i, j)
------------------ ORIG -------------------
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) + 2.0d0 * vbar(i,j)
    dvom(i,j) = dvom(i,j) + 0.5d0 * vrhs(i,j)
  enddo
enddo

-------------------------------------------
'''

##########################################################
KERNEL_4='''================= KERNEL ==================
IN    : cff1, cff2, cff3, drhs, dvom, om_v, vbar, vrhs
OUT   : dvom, vrhs
MASKS : (i, j - 1), (i, j)
------------------ ORIG -------------------
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
    dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo

-------------------------------------------

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j)
    vrhs(i,j) = vrhs(i,j) + cff2 * vbar(i,j)
    vrhs(i,j) = vrhs(i,j) + cff3 * vbar(i,j)
    dvom(i,j) = 0.5d0
    dvom(i,j) = dvom(i,j) * (drhs(i,j) + drhs(i,j - 1))
    dvom(i,j) = dvom(i,j) * om_v(i,j)
    dvom(i,j) = dvom(i,j) * vrhs(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------

------------ BEFORE FUSED AGAIN -----------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j)
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) + cff2 * vbar(i,j)
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) + cff3 * vbar(i,j)
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = 0.5d0
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * (drhs(i,j) + drhs(i,j - 1))
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * om_v(i,j)
  enddo
enddo
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * vrhs(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
'''

##########################################################
KERNEL_5='''================= KERNEL ==================
IN    : cff1
OUT   : vrhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : vbar, vrhs
OUT   : vrhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) * vbar(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : cff2, vbar, vrhs
OUT   : vrhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) + cff2 * vbar(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : cff3, vbar, vrhs
OUT   : vrhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = vrhs(i,j) + cff3 * vbar(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : 
OUT   : dvom
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = 0.5d0
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : drhs, dvom
OUT   : dvom
MASKS : (i, j - 1), (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * (drhs(i,j) + drhs(i,j - 1))
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : dvom, om_v
OUT   : dvom
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * om_v(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : dvom, vrhs
OUT   : dvom
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = dvom(i,j) * vrhs(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
'''

##########################################################
GRAPH_1='''digraph CROCO {
   cff1 -> vrhs
   vbar -> vrhs
   cff2 -> vrhs
   cff3 -> vrhs
   drhs -> dvom
   om_v -> dvom
   vrhs -> dvom
}
'''

##########################################################
GRAPH_2='''digraph CROCO
{
node_0 [label="vrhs", color="red"]
node_1 [label="cff1", color="green"]
node_1 -> node_0 
node_2 [label="vbar", color="green"]
node_2 -> node_0 
node_3 [label="cff2", color="green"]
node_3 -> node_0 
node_2 -> node_0 
node_4 [label="cff3", color="green"]
node_4 -> node_0 
node_2 -> node_0 
node_5 [label="dvom", color="red"]
node_6 [label="drhs", color="green"]
node_6 -> node_5 
node_6 -> node_5 
node_7 [label="om_v", color="green"]
node_7 -> node_5 
node_0 -> node_5 
}
'''

##########################################################
CODE_3='''
subroutine test1(cff1, cff2, cff3, vbar, drhs)
    REAL      :: vbar(0:10, 0:20)
    REAL      :: drhs(0:10, 0:20)
    REAL      :: cff1
    REAL      :: cff2
    REAL      :: cff3
    INTEGER*4 :: i
    INTEGER*4 :: j

    do j = 0, 20, 1
        do i = 0, 10, 1
            drhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
        enddo
    enddo
end

subroutine test2(cff1, cff2, cff3, vrhs, vbar, dvom, drhs, om_v)
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
        enddo
    enddo

    CALL test1(cff1, cff2, cff3, vrhs, drhs)

    do j = 0, 20, 1
        do i = 0, 10, 1
            dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
        enddo
    enddo
end
'''

##########################################################
CODE_3_SUMMARY='''================= KERNEL ==================
IN    : cff1, cff2, cff3, vbar
OUT   : drhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(1)
do j = 0, 20, 1
  do i = 0, 10, 1
    drhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : cff1, cff2, cff3, vbar
OUT   : vrhs
MASKS : (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(0)
do j = 0, 20, 1
  do i = 0, 10, 1
    vrhs(i,j) = cff1 * vbar(i,j) + cff2 * vbar(i,j) + cff3 * vbar(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
================= KERNEL ==================
IN    : drhs, om_v, vrhs
OUT   : dvom
MASKS : (i, j - 1), (i, j)

------------------ CODE -------------------
!$acc kernels default(present) async(2) wait(0, 1)
do j = 0, 20, 1
  do i = 0, 10, 1
    dvom(i,j) = 0.5d0 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo
!$acc end kernels

-------------------------------------------
'''

##########################################################
def test_basic_no_walk():
    kernel = extract_kernel_no_walk(CODE, free_form=True)
    assert kernel.render_summary(actual=False, fused=False) == KERNEL_1

##########################################################
def test_basic_walk():
    kernel = extract_kernel_walk(CODE, free_form=True)
    assert kernel.render_summary(actual=False, fused=False) == KERNEL_2

##########################################################
def test_two_kernels_walk():
    kernels = extract_kernels(CODE_2, free_form=True)
    assert kernels.render_summary(actual=False, fused=False) == KERNEL_3

##########################################################
def test_render_dep_graph():
    kernels = extract_kernels(CODE, free_form=True)
    assert kernels.render_dep_graph() == GRAPH_1

##########################################################
def test_render_step_graph():
    kernels = extract_kernels(CODE, free_form=True)
    assert kernels.render_step_graph(single_link=False) == GRAPH_2

##########################################################
def test_split_fuse():
    kernels = extract_kernels(CODE, free_form=True)
    kernels.make_acc_tranformation()
    kernels.make_loop_splitting()
    assert kernels.render_summary(actual=True, fused=True) == KERNEL_4

##########################################################
def test_split_acc_async():
    # parse
    psyir_tree = get_psyir_from_code(CODE, free_form=True)

    # splits
    kernels = extract_kernels_from_psyir(psyir_tree)
    kernels.make_loop_splitting()

    # second pass
    kernels = extract_kernels_from_psyir(psyir_tree)
    kernels.make_loop_splitting()

    # extract kernel again
    kernels = extract_kernels_from_psyir(psyir_tree)
    kernels.make_acc_tranformation(async_deps=True)
    
    assert kernels.render_summary(orig = False, actual=True, fused=False) == KERNEL_5

##########################################################
def test_split_acc_async_multi_func_1():
    # preparse
    sched = KernelPreScheduler()
    sched.add_code("main.py", CODE_3, free_form=True)
    sched.calculate2()
    sched.print()

    # parse
    psyir_tree = get_psyir_from_code(CODE_3, free_form=True)

    # apply
    sched.add_psyir("main.py", psyir_tree)
    sched.apply()

    # extract
    kernels = extract_kernels_from_psyir(psyir_tree)

    # calc stream without input
    kernels.calc_acc_async_wait(True)
    
    # apply acc rules
    kernels.make_acc_tranformation(async_deps=True)

    # gen
    assert kernels.render_summary(actual=True, orig=False, fused=False) == CODE_3_SUMMARY

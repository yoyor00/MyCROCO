##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import json
import copy
# psyclone
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.nodes.loop import Loop
from psyclone.transformations import ACCEnterDataTrans, ACCParallelTrans, \
                                     ACCLoopTrans, ACCRoutineTrans
from psyclone.psyir.nodes import ACCKernelsDirective, Assignment, Reference, ArrayReference, Routine, Call
from psyclone.psyir.nodes.schedule import Schedule
# internals
from ..base.types import AccessMode
from ..base.psyir_helpers import is_using_var, get_previous_2
from ..base.render_graph import RenderGraph
from ..transformations.make_loop_single_assignment import MakeLoopSingleAssignmentTrans
from ..transformations.make_loop_single_top_level_operand import MakeLoopSingleTopLevelOperandTrans
from ..transformations.recursive_loop_fusing import RecusiveLoopFusing

##########################################################
class PsycloneACCKernelsDirective(ACCKernelsDirective):
    '''
    Class representing the !$ACC KERNELS directive in the PSyIR.

    :param children: the PSyIR nodes to be enclosed in the Kernels region \
                     and which are therefore children of this node.
    :type children: List[:py:class:`psyclone.psyir.nodes.Node`]
    :param parent: the parent of this node in the PSyIR.
    :type parent: sub-class of :py:class:`psyclone.psyir.nodes.Node`
    :param bool default_present: whether or not to add the "default(present)" \
                                 clause to the kernels directive.

    :raises NotImplementedError: if default_present is False.

    '''
    def __init__(self, children=None, parent=None, default_present=True, default_async=None, default_wait=None, default_if=None):
        super().__init__(children=children, parent=parent)
        self._default_present = default_present
        self._default_async = default_async
        self._default_wait = default_wait
        self._default_if = default_if

    def __eq__(self, other):
        '''
        Checks whether two nodes are equal. Two ACCKernelsDirective nodes are
        equal if their default_present members are equal.

        :param object other: the object to check equality to.

        :returns: whether other is equal to self.
        :rtype: bool
        '''
        is_eq = super().__eq__(other)
        is_eq = is_eq and self.default_present == other.default_present and self._default_async == other._default_async
        is_eq = is_eq and self._default_wait == other._default_wait
        is_eq = is_eq and self._default_if == other._default_if

        return is_eq
   

    def begin_string(self):
        '''Returns the beginning statement of this directive, i.e.
        "acc kernels ...". The backend is responsible for adding the
        correct directive beginning (e.g. "!$").

        :returns: the beginning statement for this directive.
        :rtype: str

        '''
        result = "acc kernels"
        if self._default_present:
            result += " default(present)"
        if self._default_async != None:
            result += f" async({self._default_async})"
        if self._default_wait != None:
            value = ', '.join(str(value) for value in set(self._default_wait))
            result += f" wait({value})"
        if self._default_if != None:
            result += f" if({self._default_if})"

        return result

##########################################################
class Kernel:
    nextKernelId = 0

    def __init__(self, root_node: Loop):
        self.root_node = root_node
        self.orig = self.root_node.copy()
        self.before_fuse = self.orig
        self.out_vars = []
        self.in_vars = []
        self.masks = []
        self.masks_nodes = []
        self.cur_var = None
        self.var_deps = {}
        self.loop_vars = []
        self.name = f"kernel-{Kernel.nextKernelId}"
        Kernel.nextKernelId += 1
        root_node.poseidon_kernel = self

        self.acc_async_stream = 1
        self.acc_async_wait = None

        self.extract_kernel_infos()

    def make_acc_kernel(self, make_acc_in_kernel = True, stream=1, wait=None, if_clause=None):
        parent = self.root_node.parent
        start_index = self.root_node.position

        # Create a directive containing the nodes in node_list and insert it.
        directive = PsycloneACCKernelsDirective(
            parent=self.root_node.parent, children=[self.root_node.detach()],
            default_present=True, default_async=stream, default_wait=wait, default_if=if_clause)

        parent.children.insert(start_index, directive)
        if make_acc_in_kernel:
            self.root_node = directive

        #ACCKernelsTrans().apply(self.root_node, {"default_present": True})

    def make_loop_splitting(self):
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print(FortranWriter()(self.root_node))
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        self.make_loop_splitting_assignements()

        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print(FortranWriter()(self.root_node))
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        self.make_loop_splitting_first_level_operands()

        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print(FortranWriter()(self.root_node))
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        self.make_loop_splitting_assignements()

        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print(FortranWriter()(self.root_node))
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        self.fuse_again()

    def fuse_again(self):
        self.before_fuse = self.root_node.copy()
        loops = self.root_node.walk(Loop, stop_type=Loop)
        RecusiveLoopFusing().apply(loops)

    def make_loop_splitting_assignements(self):
        #assert isinstance(self.root_node, RegionDirective)
        trans = MakeLoopSingleAssignmentTrans()
        loops = self.root_node.walk(Loop, stop_type=Loop)
        for loop in loops:
            trans.apply(loop)

    def make_loop_splitting_first_level_operands(self):
        trans = MakeLoopSingleTopLevelOperandTrans()
        loops = self.root_node.walk(Loop, stop_type=Loop)
        for loop in loops:
            trans.apply(loop)

    def push_loop_var(self, name:str):
        if (not name in self.loop_vars):
            self.loop_vars.append(name)

    def extract_kernel_infos(self):
        # extract the loop vars
        loop: Loop
        for loop in self.root_node.walk(Loop):
            self.push_loop_var(loop.variable.name)

        assign: Assignment
        # handle assignments
        for assign in self.root_node.walk(Assignment):
            self.push_var(assign.lhs.name, AccessMode.WRITE)
            for var in assign.rhs.walk(Reference):
                if var.name in self.loop_vars:
                    pass
                else:
                    self.push_var(var.name, AccessMode.READ)
            self.cur_var = None

        # handle references        
        self.cur_var = None
        var: Reference
        for var in self.root_node.walk(Reference):
            # skip assign write & loop vars
            
            if isinstance(var.parent, Assignment) and var.parent.position == 0:
                pass
            elif var.name in self.loop_vars:
                pass
            else:
                self.push_var(var.name, AccessMode.READ)

        # extract stencils
        array: ArrayReference
        for array in self.root_node.walk(ArrayReference):
            if is_using_var(array.children[0], 'i'):
                stencil_i = array.children[0]
                if len(array.children) >= 2 and is_using_var(array.children[1], 'j'):
                    stencil_j = array.children[1]
                    mask = (stencil_i, stencil_j)
                    if not mask in self.masks_nodes:
                        self.masks_nodes.append(mask)
    
    def push_var(self, name: str, mode: AccessMode):
        assert mode != AccessMode.UNDEFINED
        if mode == AccessMode.WRITE:
            self.out_vars.append(name)
            self.cur_var = name
            self.var_deps[self.cur_var] = []
        elif mode == AccessMode.READ:
            self.in_vars.append(name)
            if self.cur_var != None:
                self.var_deps[self.cur_var].append(name)

    def analyse_stencil_mask(self, array_node, mode):
        self.push_var(array_node.name, mode)
        if is_using_var(array_node.children[0], 'i'):
            stencil_i = array_node.children[0]
            if is_using_var(array_node.children[1], 'j'):
                stencil_j = array_node.children[1]
                mask = (stencil_i, stencil_j)
                if not mask in self.masks_nodes:
                    self.masks_nodes.append(mask)

    @staticmethod
    def make_uniq_sorted(lst: list) -> list:
        res = list(set(lst))
        res.sort()
        return res

    @staticmethod
    def mask_to_string(mask: tuple):
        writer = FortranWriter()
        indice_0 = writer(mask[0])
        indice_1 = writer(mask[1])
        return f"({indice_0}, {indice_1})"

    def render_summary(self, orig=True, actual=True, fused=True) -> str:
        writer = FortranWriter()

        # convert masks
        seen = []
        for mask in self.masks_nodes:
            if not mask in seen:
                self.masks.append(self.mask_to_string(mask))
                seen.append(mask)

        out_vars = ', '.join(self.make_uniq_sorted(self.out_vars))
        in_vars = ', '.join(self.make_uniq_sorted(self.in_vars))
        masks = ', '.join(self.make_uniq_sorted(self.masks))
        if actual:
            code = writer(self.root_node)
        if orig:
            code_orig = writer(self.orig)
        if fused:
            code_fused = writer(self.before_fuse)

        out = f'''================= KERNEL ==================
IN    : {in_vars}
OUT   : {out_vars}
MASKS : {masks}
'''
        if orig:
            out += f'''------------------ ORIG -------------------
{code_orig}
-------------------------------------------
'''
        if actual:
            out += f'''
------------------ CODE -------------------
{code}
-------------------------------------------
'''
        if fused:
            out += f'''
------------ BEFORE FUSED AGAIN -----------
{code_fused}
-------------------------------------------
'''

        return out

    def apply_open_acc_try(self):
        #print(self.root_node.view())
        loop_trans = ACCLoopTrans()
        parallel_trans = ACCParallelTrans()
        enter_data_trans = ACCEnterDataTrans()
        rtrans = ACCRoutineTrans()

        schedule = self.root_node.loop_body
        for loop in schedule.loops():
            loop_trans.apply(loop)
            # The loop is now the child of the Directive's Schedule
            parallel_trans.apply(loop.parent.parent)
        enter_data_trans.apply(schedule)        

    def schedule(self, async_deps, var_stream, var_last_stream, streams):
        # find aync stream to be used
        selected = None
        wait = []
        if async_deps:
            #selected = 0
            for var in self.out_vars:
                if var in var_stream:
                    selected = var_stream[var]
                    break
            for var in self.in_vars + self.out_vars:
                if var in var_last_stream:
                    if var == 'ubar':
                        print(f"Wait {var} -> {var_last_stream[var]}")
                    wait.append(var_last_stream[var])
                else:
                    pass
                    if var == 'ubar':
                        print(f"Wait unknown {var}")

            if selected == None:
                selected = len(streams)
                streams.append(selected)
                wait.append(selected)
            
            for var in self.out_vars:
                if var == 'ubar':
                    print(f"{var} => {selected}")
                var_last_stream[var] = selected
                #print(f"Set {var} -> {selected}")

            #print(wait)
        else:
            selected = 1
        
        # compute final wait
        if len(wait) == 0:
            wait = None
        else:
            wait = set(wait)
            lst = list(wait)
            lst.sort()
            wait = lst
            #selected = lst[0]
            if selected in wait:
                wait.remove(selected)
            if len(wait) == 0:
                wait = None

        # store
        self.acc_async_stream = selected
        self.acc_async_wait = wait

        # save
        self.root_node.poseidon_end_state = {
            'var_last_stream': copy.deepcopy(var_last_stream),
            'streams': copy.deepcopy(streams),
            'var_stream': copy.deepcopy(var_stream),
        }

##########################################################
class KernelList:
    def __init__(self):
        self.kernels = []

    def push_kernel(self, kernel: Kernel):
        self.kernels.append(kernel)

    def get(self, id: int) -> Kernel:
        return self.kernels[id]
    
    def apply_sched_acc_async_wait(self):
        for kernel in self.kernels:
            loop = kernel.root_node.walk(Loop, stop_type=Loop)[0]
            kernel.acc_async_stream = loop.poseidon_acc_async_stream
            kernel.acc_async_wait = loop.poseidon_acc_async_wait

    def calc_acc_async_wait(self, async_deps, existing_streams=[]):
        streams = existing_streams
        var_last_stream = {}
        current_routine = ''
        for kernel in self.kernels:
            selected = None
            wait = []

            # apply func start mode
            routine = kernel.root_node.ancestor(Routine)
            if routine.name != current_routine:
                current_routine = routine.name
                try:
                    var_stream = routine.poseidon_previous_state['var_stream']
                    streams = routine.poseidon_previous_state['streams']
                    var_last_stream = routine.poseidon_previous_state['var_last_stream']
                    print(f">>>>>>>>> routine={current_routine}")
                    print(streams)
                    print(var_last_stream)
                except:
                    streams = []
                    var_last_stream = {}
                    var_stream = {}

            # is preceded by al function call
            call = get_previous_2(kernel.root_node, (Loop, Call))
            if call != None:
                try:
                    streams = call.poseidon_end_state['streams']
                    var_last_stream = call.poseidon_end_state['var_last_stream']
                    var_stream = call.poseidon_end_state['var_stream']
                    print(f">>>>>>>>> call={call.name}")
                    print(streams)
                    print(var_last_stream)
                except:
                    pass

            # schdule in kernel
            kernel.schedule(async_deps, var_stream, var_last_stream, streams)

            # save for pre analysis
            routine.poseidon_end_state = kernel.root_node.poseidon_end_state
        
        # ret
        return streams

    def make_acc_tranformation(self, make_acc_in_kernel = True, async_deps = False, if_clause = None):
        #self.calc_acc_async_wait(async_deps)
        for kernel in self.kernels:
            kernel.make_acc_kernel(make_acc_in_kernel, stream=kernel.acc_async_stream, wait=kernel.acc_async_wait, if_clause = if_clause)

    def merge_joinable_kernels(self):
        first_acc_to_merge = None
        first_kernel_to_merge = None
        for kernel in self.kernels:
            if first_acc_to_merge == None:
                assert isinstance(kernel.root_node, Loop)
                assert isinstance(kernel.root_node.parent.parent, ACCKernelsDirective)
                first_acc_to_merge = kernel.root_node.parent.parent
                first_kernel_to_merge = kernel
            elif first_acc_to_merge.parent == kernel.root_node.parent.parent.parent:
                #print(first_acc_to_merge.view())
                schedule = kernel.root_node.parent
                acc = schedule.parent
                assert isinstance(kernel.root_node, Loop)
                assert isinstance(schedule, Schedule)
                assert isinstance(first_acc_to_merge, ACCKernelsDirective)
                assert isinstance(acc, ACCKernelsDirective)
                print(f"merge {first_kernel_to_merge.name} and {kernel.name}")
                print(f"merge {first_acc_to_merge.position} and {acc.position}")
                #print("----------------------------")
                #print(first_acc_to_merge.view())
                #print("----------------------------")
                #print(acc.view())
                print("----------------------------")
                if first_acc_to_merge.position == acc.position - 1:
                    schedule.detach()
                    childs = []
                    for child in schedule.children:
                        child.detach()
                        childs.append(child)
                        first_acc_to_merge.children[0].children.append(child)
                    acc.detach()
                else:
                    first_acc_to_merge = acc
                    first_kernel_to_merge = kernel
            else:
                first_acc_to_merge = kernel.root_node.parent.parent
                first_kernel_to_merge = kernel
                

    def make_loop_splitting(self):
        for kernel in self.kernels:
            kernel.make_loop_splitting()

    def make_loop_fusing(self):
        for kernel in self.kernels:
            kernel.fuse_again()

    def render_summary(self, orig=True, actual=True, fused=True):
        out = ""
        for kernel in self.kernels:
            out += kernel.render_summary(orig=orig, actual=actual, fused=fused)
        return out

    def render_dep_graph(self, single_link = True) -> str:
        out = "digraph CROCO {\n"
        seen = []
        for kernel in self.kernels:
            for out_var, in_vars in kernel.var_deps.items():
                for in_var in in_vars:
                    tmp = f"{in_var} -> {out_var}"
                    if not tmp in seen:
                        out += f"   {tmp}\n"
                        if single_link:
                            seen.append(tmp)
        return out + "}\n"
    
    def render_dep_graph_vars(self, single_link = True) -> str:
        graph = RenderGraph(digraph=True)
        seen = []
        vars = {}
        var_versions = {}
        for kernel in self.kernels:
            out_vars = []
            for out_var in kernel.var_deps:
                out_vars.append(out_var)
            for out_var, in_vars in kernel.var_deps.items():
                version=1
                if out_var in var_versions:
                    version = var_versions[out_var] + 1
                new_var = graph.add_node(f"{out_var} (v{version})", None)
                for in_var in in_vars + out_vars:
                    if in_var in vars:
                        in_kernel = vars[in_var]
                        tmp = f"{in_kernel.name} -> {new_var.name}"
                        if not tmp in seen:
                            if out_var == in_var:
                                graph.add_link(in_kernel, new_var, style='style="dashed"')
                            else:
                                graph.add_link(in_kernel, new_var)
                            if single_link:
                                seen.append(tmp)
                vars[out_var] = new_var
                var_versions[out_var] = version
                
        return graph.render()
    
    def render_dep_graph_kernels(self, single_link = True, exec_order_links = False) -> str:
        graph = RenderGraph(digraph=True)
        seen = []
        vars = {}
        prev_kernel_node = None
        for id, kernel in enumerate(self.kernels):
            lst_in = []
            lst_out = []
            for out_var, in_vars in kernel.var_deps.items():
                if not out_var in lst_out:
                    lst_out.append(out_var)
                for in_var in in_vars:
                    if not in_var in lst_in:
                        lst_in.append(in_var)
            out_str = ','.join(set(lst_out))
            in_str = ','.join(set(lst_in))

            kernel_node = graph.add_kernel(f"Kernel {id}\nIN: {in_str}\nOUT: {out_str}")
            var_chain = kernel_node
            out_vars = []
            for out_var in kernel.var_deps:
                out_vars.append(out_var)
            for out_var, in_vars in kernel.var_deps.items():
                for in_var in in_vars + out_vars:
                    if in_var in vars:
                        in_kernel = vars[in_var]
                        tmp = f"{in_kernel.name} -> {kernel_node.name}"
                        if not tmp in seen:
                            graph.add_link(in_kernel, kernel_node)
                            if single_link:
                                seen.append(tmp)
            for var in set(lst_out):
                vars[var] = kernel_node

            if id > 0 and exec_order_links:
                graph.add_link(prev_kernel_node, kernel_node, 'color=red')
            prev_kernel_node = kernel_node

        return graph.render()

    def render_step_graph(self, single_link = False) -> str:
        seen = []
        graph = RenderGraph(digraph=True)
        last_node = {}
        for kernel in self.kernels:
            if len(kernel.loop_vars) == 2:
                for out_var, in_vars in kernel.var_deps.items():
                    out_node = graph.add_node(out_var, mode=AccessMode.WRITE)
                    for in_var in in_vars:
                        if in_var in last_node:
                            in_node = last_node[in_var]
                        else:
                            in_node = graph.add_node(in_var, mode=AccessMode.READ)
                            last_node[in_var] = in_node
                        tmp = f"{out_node} -> {in_node}"
                        if not tmp in seen:
                            graph.add_link(in_node, out_node)
                            if single_link:
                                seen.append(tmp)
                    last_node[out_var] = out_node
        return graph.render()

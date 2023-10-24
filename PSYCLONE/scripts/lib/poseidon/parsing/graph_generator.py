##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
#PSyIR internal representation
from psyclone.psyir.nodes import Routine, Reference, Assignment, ArrayReference, \
     Loop, Call, IfBlock, Schedule, CodeBlock, Node
# internal
from .walker import WalkerCallbackInterface, WalkerCutRecurse
from ..base.types import AccessMode
from ..base.render_graph import RenderGraph
from ..base.psyir_helpers import extract_var_ref_list

##########################################################
class GraphGenerator(WalkerCallbackInterface):
    def __init__(self) -> None:
        self.stack = []
        self.graph = RenderGraph()

    def start_walking(self):
        pass

    def on_routine(self, routine: Routine, enter: bool) -> None:
        if enter:
            self.graph.open_subgraph(routine.name)
            graph_node = self.graph.add_node(f"ROOT {routine.name}")
            self.stack.append(graph_node)
        else:
            self.graph.close_subgraph()
            self.stack.pop()

    def on_loop(self, loop: Loop, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_loop(loop.variable.name, parent_name=self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_loop_limits(self, loop: Loop, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("limits", parent_name=self.stack[-1])
            self.stack.append(graph_node)

            deps = []
            deps += extract_var_ref_list(loop.start_expr,0)
            deps += extract_var_ref_list(loop.stop_expr,0)
            deps += extract_var_ref_list(loop.step_expr,0)

            for limit in deps:
                self.graph.add_node(limit['var'] + " (" + limit['mode'] + ")", parent_name=self.stack[-1], mode=limit['mode'])

            raise WalkerCutRecurse
        else:
            self.stack.pop()

    def on_loop_body(self, loop: Loop, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node('body', parent_name = self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_array_reference(self, array: ArrayReference, access_mode: AccessMode, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("Array: " + array.name + f" ({access_mode})", parent_name=self.stack[-1], mode=access_mode)
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_reference(self, variable: Reference, access_mode: AccessMode, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node(variable.name + f" ({access_mode})", parent_name=self.stack[-1], mode=access_mode)
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_call(self, call: Call, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_call(call.routine, parent_name=self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_schedule(self, schedule: Schedule, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_schedule(parent_name=self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_if(self, if_block: IfBlock, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_if(parent_name=self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_if_cond(self, if_block: IfBlock, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("condition", parent_name = self.stack[-1])
            self.stack.append(graph_node)

            deps = extract_var_ref_list(if_block.condition, 0)
            for var in deps:
                self.graph.add_node(var['var'] + " (" + var['mode'] + ")", parent_name=graph_node, mode=var['mode'])

            raise WalkerCutRecurse
        else:
            self.stack.pop()

    def on_if_body(self, if_block: IfBlock, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("body", parent_name = self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_if_else(self, if_block: IfBlock, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("else", parent_name = self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_assignment(self, assign: Assignment, enter: bool) -> None:
        if enter:
            graph_node = self.graph.add_node("=", parent_name = self.stack[-1])
            self.stack.append(graph_node)
        else:
            self.stack.pop()

    def on_code_block(self, block: CodeBlock, enter: bool) -> None:
        pass

    def on_others(self, node: Node, enter: bool) -> None:
        pass

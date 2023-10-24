##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
#PSyIR internal representation
from psyclone.psyir.nodes import Routine, Reference, Assignment, ArrayReference, \
    FileContainer, BinaryOperation, Literal, Return, IfBlock, Schedule, \
    UnaryOperation, CodeBlock, Node, Loop, Call
from psyclone.nemo import NemoKern
# internal
from ..base.types import AccessMode

##########################################################
class WalkerInterupt(Exception):
    pass

##########################################################
class WalkerCutRecurse(Exception):
    pass

##########################################################
class WalkerCallbackInterface:
    def __init__(self) -> None:
        pass

    def start_walking(self):
        pass

    def on_routine(self, routine: Routine, enter: bool) -> None:
        pass

    def on_loop(self, loop: Loop, enter: bool) -> None:
        pass

    def on_loop_limits(self, loop: Loop, enter: bool) -> None:
        pass

    def on_loop_body(self, loop: Loop, enter: bool) -> None:
        pass

    def on_array_reference(self, array: ArrayReference, access_mode: AccessMode, enter: bool) -> None:
        pass

    def on_reference(self, variable: Reference, access_mode: AccessMode, enter: bool) -> None:
        pass

    def on_call(self, call: Call, enter: bool) -> None:
        pass

    def on_schedule(self, schedule: Schedule, enter: bool) -> None:
        pass

    def on_if(self, if_block: IfBlock, enter: bool) -> None:
        pass

    def on_if_cond(self, if_block: IfBlock, enter: bool) -> None:
        pass

    def on_if_body(self, if_block: IfBlock, enter: bool) -> None:
        pass

    def on_if_else(self, if_block: IfBlock, enter: bool) -> None:
        pass

    def on_assignment(self, assign: Assignment, enter: bool) -> None:
        pass

    def on_code_block(self, block: CodeBlock, enter: bool) -> None:
        pass

    def on_others(self, node: Node, enter: bool) -> None:
        pass

##########################################################
class Walker:
    def __init__(self, callback: WalkerCallbackInterface) -> None:
        self.callback = callback

    def walk(self, psyir_tree: Node) -> None:
        self.callback.start_walking()
        try:
            self.walk_on_node(psyir_tree)
        except WalkerInterupt:
            pass

    def walk_on_childs(self, node: Node) -> None:
        for id, child in enumerate(node.children):
            self.walk_on_node(child, id)


    @staticmethod
    def get_access_mode(node: Node, id: int) -> AccessMode:
        if isinstance(node.parent, Assignment) and id == 0:
            return AccessMode.WRITE
        else:
            return AccessMode.READ

    def walk_on_node(self, node: Node, id: int = 0) -> None:
        try:
            if isinstance(node, Routine):
                self.callback.on_routine(node, True)
                self.walk_on_childs(node)
                self.callback.on_routine(node, False)
            elif isinstance(node, Loop):
                self.callback.on_loop(node, True)
                
                try:
                    self.callback.on_loop_limits(node, True)
                    self.walk_on_node(node.start_expr, 0)
                    self.walk_on_node(node.stop_expr, 0)
                    self.walk_on_node(node.step_expr, 0)
                    self.callback.on_loop_limits(node, False)
                except WalkerCutRecurse:
                    pass

                try:
                    self.callback.on_loop_body(node, True)
                    self.walk_on_node(node.loop_body, 0)
                    self.callback.on_loop_body(node, False)
                except WalkerCutRecurse:
                    pass           
                
                self.callback.on_loop(node, False)
            elif isinstance(node, ArrayReference):
                mode = self.get_access_mode(node, id)
                self.callback.on_array_reference(node, mode, True)
                self.walk_on_childs(node)
                self.callback.on_array_reference(node, mode, False)
            elif isinstance(node, Reference):
                if node.name == 'i' or node.name == 'j':
                    pass
                mode = self.get_access_mode(node, id)
                self.callback.on_reference(node, mode, True)
                self.walk_on_childs(node)
                self.callback.on_reference(node, mode, False)
            elif isinstance(node, Call):
                self.callback.on_call(node, True)
                self.walk_on_childs(node)
                self.callback.on_call(node, False)
            elif isinstance(node, Schedule):
                self.callback.on_schedule(node, True)
                self.walk_on_childs(node)
                self.callback.on_schedule(node, False)
            elif isinstance(node, IfBlock):
                self.callback.on_if(node, True)

                try:
                    self.callback.on_if_cond(node, True)
                    self.walk_on_node(node.condition, 0)
                    self.callback.on_if_cond(node, False)
                except WalkerCutRecurse:
                    pass

                try:
                    self.callback.on_if_body(node, True)
                    self.walk_on_node(node.if_body, 1)
                    self.callback.on_if_body(node, False)
                except WalkerCutRecurse:
                    pass

                try:
                    if 2 in node.children:
                        self.callback.on_if_else(node, True)
                        self.walk_on_node(node.else_body, 2)
                        self.callback.on_if_else(node, False)
                except WalkerCutRecurse:
                    pass
                
                self.callback.on_if(node, False)
            elif isinstance(node, Assignment):
                self.callback.on_assignment(node, True)
                self.walk_on_childs(node)
                self.callback.on_assignment(node, False)
            elif isinstance(node, CodeBlock):
                self.callback.on_code_block(node, True)
                self.walk_on_childs(node)
                self.callback.on_code_block(node, False)
            elif isinstance(node, FileContainer) or isinstance(node, BinaryOperation) or isinstance(node, Literal) or isinstance(node, Return) or isinstance(node, UnaryOperation):
                self.callback.on_others(node, True)
                self.walk_on_childs(node)
                self.callback.on_others(node, False)
            elif isinstance(node, NemoKern):
                self.walk_on_childs(node)
            else:
                typename = type(node)
                raise Exception(f"Unsupported type : {typename} !")
        except WalkerCutRecurse:
            pass

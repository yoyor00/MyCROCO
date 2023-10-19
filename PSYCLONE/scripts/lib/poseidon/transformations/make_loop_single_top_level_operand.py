# psyir
from psyclone.psyir.nodes.loop import Loop

# transf
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes.assignment import Assignment
from psyclone.psyir.nodes.array_reference import ArrayReference
from psyclone.psyir.nodes.reference import Reference
from psyclone.psyir.nodes.schedule import Schedule
from psyclone.psyGen import Transformation
from psyclone.psyir.nodes.operation import BinaryOperation

class MakeLoopSingleTopLevelOperandTrans(Transformation):
    def __str__(self):
        return "If there is multiple top level + or - or * operands, we split the operation in multiple lines."

    @property
    def name(self):
        return 'MakeLoopSingleTopLevelOperandTrans'

    @staticmethod
    def _get_inner_loop(loop: Loop):
        assert isinstance(loop, (Loop, Schedule))
        inner_loop = loop
        #print(loop.view())
        while isinstance(inner_loop, (Loop, Schedule)):
            if isinstance(inner_loop, Schedule):
                if len(inner_loop.children) == 0:
                    return inner_loop
                else:
                    inner_loop = inner_loop.children[0]
            elif isinstance(inner_loop, Loop):
                inner_loop = inner_loop.loop_body
        return inner_loop.parent

    def validate(self, top_loop, options = {}):
        # top node should be a loop
        if not isinstance(top_loop, Loop):
            raise TransformationError(f"Get invalid Loop as top loop while running {self.name} : got {type(top_loop)} !")

        # we handle only direct loop childing, not calculation except child loop in intermediate loops
        for loop in top_loop.walk(Loop)[-1:]:
            if not isinstance(loop.parent.parent, Loop):
                raise TransformationError(f"Support direct chaining between loops, got parent : {loop.parent}")
            if not len(loop.parent.children) == 1:
                raise TransformationError(f"Unsupported complex loop, not getting Loop -> Loop -> Action, but habing unsupported more complex schema.")

        # we currently do not handle write to scalar values in child loops
        for scalar_var in top_loop.walk(Reference):
            # this is assignement to an array, can skip
            if isinstance(scalar_var, ArrayReference):
                continue

            # if parent is assigned & we are lhs (assigne to, so position 0)
            scalar_parent = scalar_var.parent
            if isinstance(scalar_parent, Assignment) and scalar_var.position == 0:
                raise TransformationError(f"We currently do not handle splitting loops with scalar assignement in it ({scalar_var.name}).")

        # body line is an assignment
        inner_loop = self._get_inner_loop(top_loop)
        if len(inner_loop.children) != 1:
            raise TransformationError("Have multiple child operations, need to have only one !")
        
        # check op type
        op = inner_loop.children[0]
        if not isinstance(op, Assignment):
            raise TransformationError("Inner operation should be an assignement !")
        if not isinstance(op.children[0], ArrayReference):
            raise TransformationError("Assignement should be directed to an array reference, not a scalar !")

        # nothing to do
        if not isinstance(op.children[1], BinaryOperation):
            raise TransformationError("Assignement do not contain binary operation as root element, nothing to transform !")

    @staticmethod
    def _is_both_mul(operation1: BinaryOperation, operation2: BinaryOperation):
        return operation1.operator == BinaryOperation.Operator.MUL and operation1.operator == operation2.operator

    @staticmethod
    def _is_both_add_sub(operation1: BinaryOperation, operation2: BinaryOperation):
        if operation1.operator != BinaryOperation.Operator.ADD and operation1.operator != BinaryOperation.Operator.SUB:
            return False
        if operation2.operator != BinaryOperation.Operator.ADD and operation2.operator != BinaryOperation.Operator.SUB:
            return False
        return True

    def apply(self, top_loop, options = {}):
        # validate
        self.validate(top_loop, options=options)

        # jump to inner loop to extract content
        inner_loop = self._get_inner_loop(top_loop)

        # abort if has multi line
        if len(inner_loop.children) != 1:
            return

        # get to operation
        op = inner_loop.children[0]

        # safe checj
        assert isinstance(op, Assignment)
        assert isinstance(op.children[0], ArrayReference)

        # nothing to do
        if not isinstance(op.children[1], BinaryOperation):
            return

        # extract infos
        root_op = op.rhs
        op_cursor = root_op
        assign = op.lhs
        is_self_ref = False
        lhs = None
        rhs = None

        # loop while we have binary op we can split
        while isinstance(op_cursor, BinaryOperation) and (self._is_both_mul(root_op, op_cursor) or self._is_both_add_sub(root_op, op_cursor)):
            lhs = op_cursor.children[0]
            rhs = op_cursor.children[1]
            if rhs == assign:
                is_self_ref
            else:
                inner_loop.addchild(Assignment.create(assign.copy(), 
                                                    BinaryOperation.create(
                                                        op_cursor.operator, assign.copy(), rhs.copy())), 0)
            op_cursor = lhs
        if lhs != None and assign != lhs:
            if is_self_ref:
                inner_loop.addchild(Assignment.create(assign.copy(), 
                                                    BinaryOperation.create(
                                                        op_cursor.operator, assign.copy(), lhs.copy())), 0)
            else:
                inner_loop.addchild(Assignment.create(assign.copy(), lhs.copy()), 0)
        
        if lhs != None:
            op.detach()

    @property
    def name(self):
        '''
        :returns: the name of this class.
        :rtype: str
        '''
        return self.__class__.__name__

# psyir
from psyclone.psyir.nodes.loop import Loop

# transf
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes.assignment import Assignment
from psyclone.psyir.nodes.array_reference import ArrayReference
from psyclone.psyir.nodes.reference import Reference
from psyclone.psyir.nodes.schedule import Schedule
from psyclone.psyGen import Transformation, InlinedKern

class MakeLoopSingleAssignmentTrans(Transformation):
    def __str__(self):
        return "If there is multiple affectation lines, split then to get one in each loop set."

    @property
    def name(self):
        return 'MakeLoopSingleLineTrans'

    @staticmethod
    def _get_inner_loop(loop: Loop):
        assert isinstance(loop, (Loop, Schedule))
        inner_loop = loop
        #print("----------------------------------")
        #print(loop.view())
        while isinstance(inner_loop, (Loop, Schedule, InlinedKern)):
            if isinstance(inner_loop, (Schedule, InlinedKern)):
                if len(inner_loop.children) == 0:
                    return inner_loop
                else:
                    inner_loop = inner_loop.children[0]
            elif isinstance(inner_loop, Loop):
                inner_loop = inner_loop.loop_body
        #print(inner_loop.parent)
        #print("----------------------------------")
        return inner_loop.parent

    def validate(self, top_loop, options = {}):
        # top node should be a loop
        if not isinstance(top_loop, Loop):
            raise TransformationError(f"Get invalid Loop as top loop while running {self.name} : got {type(top_loop)} !")

        # we handle only direct loop childing, not calculation except child loop in intermediate loops
        for loop in top_loop.walk(Loop)[-1:]:
            if not isinstance(loop.parent.parent, Loop):
                raise TransformationError(f"Support only direct chaining between loops, got parent : {type(loop.parent.parent)}")
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

    def apply(self, top_loop, options = {}):
        # validate
        self.validate(top_loop, options=options)

        # get parent node in which to create the sub loops
        parent_node = top_loop.parent
        top_loop_position = top_loop.position

        # jump to inner loop to extract content
        inner_loop = self._get_inner_loop(top_loop)

        # alredy have a single child, nothing to split
        if len(inner_loop.children) == 1:
            return

        # extract content
        ops = []
        to_detach = []
        for op in inner_loop.children:
            if not isinstance(op, Assignment):
                TransformationError(f"Expect Assignement, get {type(op)}")
            assert isinstance(op, Assignment)
            assert isinstance(op.children[0], ArrayReference)
            to_detach.append(op)

        # detach all child ops to keep the nested loops as template to copy
        for op in to_detach:
            ops.append(op)
            op.detach()

        # rebuild
        new_loops = [top_loop]
        pos = top_loop_position
        for i in range(len(ops) - 1):
            new_empty_loop = top_loop.copy()
            new_loops.append(new_empty_loop)

        # fill
        for i, op in enumerate(ops):
            inner_new_loop = self._get_inner_loop(new_loops[i])
            inner_new_loop.addchild(op)

        # attach loops
        for i in range(len(new_loops) - 1):
            parent_node.addchild(new_loops[i+1], pos + i + 1)

    @property
    def name(self):
        '''
        :returns: the name of this class.
        :rtype: str
        '''
        return self.__class__.__name__

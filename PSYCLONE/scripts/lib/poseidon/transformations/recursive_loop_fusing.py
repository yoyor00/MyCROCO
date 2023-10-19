# psyir
from psyclone.psyir.nodes.loop import Loop

# transf
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes.assignment import Assignment
from psyclone.psyir.nodes.array_reference import ArrayReference
from psyclone.psyir.nodes.reference import Reference
from psyclone.psyir.nodes.schedule import Schedule
from psyclone.psyGen import Transformation
from psyclone.psyir.transformations.loop_fuse_trans import LoopFuseTrans

class RecusiveLoopFusing(Transformation):
    def __str__(self):
        return "Fuse neihtboor loops."

    @property
    def name(self):
        return 'MakeLoopSingleLineTrans'

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

    def validate(self, loops, options = {}):
        # perform some checkings
        current_pos = -1
        for loop in loops:
            # they should be loop
            if not isinstance(loop, Loop):
                raise TransformationError(f"Get invalid Loop as top loop while running {self.name} : got {type(loop)} !")
            # they should have same parent
            if loop.parent != loops[0].parent:
                raise TransformationError(f"Don't get same parent for one of the loop !")
            # check they are one after the other
            if current_pos != -1 and loop.position != current_pos + 1:
                raise TransformationError(f"The loops should be neighboors without intermediate ops !")
            # check we can merge
            if current_pos != -1:
                LoopFuseTrans().validate(loops[0], loop, options)

    def apply(self, loops, options = {}):

        # validate
        self.validate(loops, options)

        # furse level &
        for i in range(len(loops) - 1):
            trans = LoopFuseTrans()
            trans.apply(loops[0], loops[i+1])

        # extract childs
        next_level_loops = []
        for loop in loops[0].loop_body.children:
            assert isinstance(loop, Loop)
            next_level_loops.append(loop)

        # fuse childs if possible
        try:
            RecusiveLoopFusing().validate(next_level_loops)
            RecusiveLoopFusing().apply(next_level_loops)
        except:
            pass

    @property
    def name(self):
        '''
        :returns: the name of this class.
        :rtype: str
        '''
        return self.__class__.__name__

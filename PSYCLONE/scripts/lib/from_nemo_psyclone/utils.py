# imported from PSyClone NEMO scripts

''' Utilities file to parallelise Nemo code. '''

from psyclone.domain.nemo.transformations import NemoAllArrayRange2LoopTrans
from psyclone.errors import InternalError
from psyclone.psyir.nodes import Loop, Assignment, Directive, CodeBlock, Call, Return
from psyclone.psyir.transformations import HoistLoopBoundExprTrans, HoistTrans
from psyclone.transformations import TransformationError, ACCKernelsTrans


def normalise_loops(
        schedule,
        unwrap_array_ranges: bool = True,
        hoist_expressions: bool = True,
        ):
    ''' Normalise all loops in the given schedule so that they are in an
    appropriate form for the Parallelisation transformations to analyse
    them.

    :param schedule: the PSyIR Schedule to transform.
    :param unwrap_array_ranges: whether to convert ranges to explicit loops.
    :param hoist_expressions: whether to hoist bounds and loop invariant \
        statements out of the loop nest.
    '''
    if unwrap_array_ranges:
        # Convert all array implicit loops to explicit loops
        explicit_loops = NemoAllArrayRange2LoopTrans()
        for assignment in schedule.walk(Assignment):
            explicit_loops.apply(assignment)

    if hoist_expressions:
        # First hoist all possible expressions
        for loop in schedule.walk(Loop):
            HoistLoopBoundExprTrans().apply(loop)

        # Hoist all possible assignments (in reverse order so the inner loop
        # constants are hoisted all the way out if possible)
        for loop in reversed(schedule.walk(Loop)):
            for statement in list(loop.loop_body):
                try:
                    HoistTrans().apply(statement)
                except TransformationError:
                    pass


def insert_explicit_loop_parallelism(
        schedule,
        region_directive_trans=None,
        loop_directive_trans=None,
        collapse: bool = True
        ):
    ''' For each loop in the schedule that doesn't already have a Directive
    as an ancestor, attempt to insert the given region and loop directives.

    :param region_directive_trans: PSyclone transformation to insert the \
        region directive.
    :param loop_directive_trans: PSyclone transformation to use to insert the \
        loop directive.
    :param collapse: whether to attempt to insert the collapse clause to as \
        many nested loops as possible.
    '''

    # Add the parallel directives in each loop
    for loop in schedule.walk(Loop):
        if loop.ancestor(Directive):
            continue  # Skip if an outer loop is already parallelised

        try:
            loop_directive_trans.apply(loop)
            # Only add the region directive if the loop was successfully
            # parallelised.
            if region_directive_trans is not None:
                region_directive_trans.apply(loop.parent.parent)
        except TransformationError as err:
            # This loop can not be transformed, proceed to next loop
            print("Loop not parallelised because:", str(err))
            continue

        if collapse:
            # Count the number of perfectly nested loops
            num_nested_loops = 0
            next_loop = loop
            while isinstance(next_loop, Loop):
                num_nested_loops += 1
                if len(next_loop.loop_body.children) > 1:
                    break
                next_loop = next_loop.loop_body.children[0]

            if num_nested_loops > 1:
                loop.parent.parent.collapse = num_nested_loops


def valid_kernel(node):
    '''
    Whether the sub-tree that has `node` at its root is eligible to be
    enclosed within an OpenACC KERNELS directive.

    :param node: the node in the PSyIR to check.
    :type node: :py:class:`psyclone.psyir.nodes.Node`

    :returns: True if the sub-tree can be enclosed in a KERNELS region.
    :rtype: bool

    '''
    excluded_node_types = (CodeBlock, Call, Return)
    return node.walk(excluded_node_types) == []


def add_kernels(children, default_present=True, async_queue=False):
    '''
    Walks through the PSyIR inserting OpenACC KERNELS directives at as
    high a level as possible.

    :param children: list of sibling Nodes in PSyIR that are candidates for \
                     inclusion in an ACC KERNELS region.
    :type children: list of :py:class:`psyclone.psyir.nodes.Node`
    :param bool default_present: whether or not to supply the \
                          DEFAULT(PRESENT) clause to ACC KERNELS directives.

    '''
    if not children:
        return

    node_list = []
    for child in children[:]:
        # Can this node be included in a kernels region?
        if not valid_kernel(child):
            try_kernels_trans(node_list, default_present, async_queue)
            node_list = []
            # It can't so go down a level and try again
            # add_kernels(child.ACCLoopDirective must have an ACCParallelDirective or ACCKechildren)
        else:
            node_list.append(child)

    try_kernels_trans(node_list, default_present, async_queue)


def try_kernels_trans(nodes, default_present, async_queue):
    '''
    Attempt to enclose the supplied list of nodes within a kernels
    region. If the transformation fails then the error message is
    reported but execution continues.

    :param nodes: list of Nodes to enclose within a Kernels region.
    :type nodes: list of :py:class:`psyclone.psyir.nodes.Node`
    :param bool default_present: whether or not to supply the \
                          DEFAULT(PRESENT) clause to ACC KERNELS directives.

    '''
    have_one = False
    for node in nodes:
        if isinstance(node, Loop) or not len(node.walk(Loop)) == 0:
            have_one = True
    if not have_one:
        return
    try:
        ACCKernelsTrans().apply(nodes, {"default_present": default_present, "async_queue": async_queue})
    except (TransformationError, InternalError) as err:
        print(f"Failed to transform nodes: {nodes}")
        print(f"Error was: {err}")

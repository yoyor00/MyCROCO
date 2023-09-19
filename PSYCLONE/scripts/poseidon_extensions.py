# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2018-2022, Science and Technology Facilities Council.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------
# Authors: R. W. Ford, A. R. Porter and S. Siso, STFC Daresbury Lab

'''A transformation script that seeks to apply OpenACC KERNELS directives to
NEMO style code. In order to use it you must first install PSyclone. See
README.md in the top-level directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s <this_script> <target_source_file>

The transformation script attempts to insert Kernels directives at the
highest possible location(s) in the schedule tree (i.e. to enclose as
much code as possible in each Kernels region).

'''

from psyclone.core import Signature
from psyclone.psyir.nodes import ACCDirective, RegionDirective, ACCStandaloneDirective
from psyclone.f2pygen import DirectiveGen
import abc
from psyclone.psyGen import Transformation
from poseidon.dsl.kernel import PsycloneACCKernelsDirective
from psyclone.psyir.nodes import ACCDataDirective, ACCDirective, \
    ACCEnterDataDirective, ACCKernelsDirective, ACCLoopDirective, \
    ACCParallelDirective, ACCRoutineDirective, Assignment, CodeBlock, \
    Directive, Loop, Node, OMPDeclareTargetDirective, \
    OMPDirective, OMPMasterDirective, \
    OMPParallelDirective, OMPParallelDoDirective, OMPSerialDirective, \
    OMPSingleDirective, OMPTaskloopDirective, PSyDataNode, Reference, \
    Return, Routine, Schedule
from utils import add_kernels, normalise_loops, \
    insert_explicit_loop_parallelism
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from psyclone.psyir.transformations import ACCUpdateTrans
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.if_block import IfBlock

class ACCWaitDirective(ACCStandaloneDirective):
    '''
    Class representing the !$ACC KERNELS directive in the PSyIR.

    :param wait_queue: Which ACC async group to wait. None to wait all.
    :type wait_queue: None or int or Signature
    '''
    def __init__(self, wait_queue=None, async_queue=None):
        # call parent
        super().__init__()

        # check
        if wait_queue != None and not isinstance(wait_queue, (int, list, Signature)):
            raise TypeError("Invalid value type as wait_group, shoule be in (None, int, Signature) !")

        # assign
        self._wait_queue = wait_queue
        self._async_queue = async_queue

    def __eq__(self, other):
        '''
        Test the equality of two directives.

        :returns: If the two directives are equals.
        :rtype: bool
        '''
        is_eq = super().__eq__(other)
        is_eq = is_eq and self._wait_queue == other._wait_queue
        is_eq = is_eq and self._async_queue == other._async_queue
        return is_eq

    @property
    def wait_queue(self):
        '''
        :returns: Define which queue to wait.
        :rtype: None or str or int
        '''
        return self._wait_queue
    
    @wait_queue.setter
    def wait_queues(self, wait_queue):
        # check
        if wait_queue != None and not isinstance(wait_queue, (int, list, Signature)):
            raise TypeError("Invalid value type as wait_group, shoule be in (None, int, Signature) !")
        
        # set
        self._wait_queue = wait_queue

    def gen_code(self, parent):
        '''
        Generate the given directive code to add it to the call tree.

        :param parent: the parent Node in the Schedule to which to add this \
                       content.
        :type parent: sub-class of :py:class:`psyclone.f2pygen.BaseGen`
        '''
        # Generate the directive
        args = ' '.join(self.begin_string().split()[2:])
        parent.add(DirectiveGen(parent, "acc", "begin", "wait", args))

    def begin_string(self):
        '''Returns the beginning statement of this directive, i.e.
        "acc wait ...". The backend is responsible for adding the
        correct directive beginning (e.g. "!$").

        :returns: the beginning statement for this directive.
        :rtype: str

        '''
        # default basic directive
        result = "acc wait"

        # wait
        wait = self._wait_queue
        if isinstance(wait, list):
            wait = ', '.join(wait)

        # handle specifying groups
        if self._wait_queue != None:
            if isinstance(self._wait_queue, Signature):
                result += f" ({self._wait_queue.var_name})"
            else:
                result += f" ({self._wait_queue})"
    
        # async
        if self._async_queue != None:
            if isinstance(self._wait_queue, Signature):
                result += f" async({self._async_queue.var_name})"
            else:
                result += f" async({self._async_queue})"

        # ok return it
        return result

class ACCSetDeviceNumDirective(ACCStandaloneDirective):
    '''
    Class representing the !$ACC KERNELS directive in the PSyIR.

    :param device_num: Which ACC async group to wait.
    :type device_num: int or Signature
    '''
    def __init__(self, device_num):
        # call parent
        super().__init__()

        # check
        if not isinstance(device_num, (int, Signature)):
            raise TypeError("Invalid value type as wait_group, shoule be in (int, Signature) !")

        # assign
        self._device_num = device_num

    def __eq__(self, other):
        '''
        Test the equality of two directives.

        :returns: If the two directives are equals.
        :rtype: bool
        '''
        is_eq = super().__eq__(other)
        is_eq = is_eq and self._device_num == other._device_num
        return is_eq

    @property
    def device_num(self):
        '''
        :returns: Define which queue to wait.
        :rtype: None or str or int
        '''
        return self._device_num
    
    @device_num.setter
    def device_num(self, device_num):
        # check
        if not isinstance(device_num, (int, Signature)):
            raise TypeError("Invalid value type as wait_group, shoule be in (None, int, Signature) !")
        
        # set
        self._device_num = device_num

    def gen_code(self, parent):
        '''
        Generate the given directive code to add it to the call tree.

        :param parent: the parent Node in the Schedule to which to add this \
                       content.
        :type parent: sub-class of :py:class:`psyclone.f2pygen.BaseGen`
        '''
        # Generate the directive
        args = ' '.join(self.begin_string().split()[2:])
        parent.add(DirectiveGen(parent, "acc", "set", "device_num", args))

    def begin_string(self):
        '''Returns the beginning statement of this directive, i.e.
        "acc wait ...". The backend is responsible for adding the
        correct directive beginning (e.g. "!$").

        :returns: the beginning statement for this directive.
        :rtype: str

        '''

        # handle specifying groups
        if isinstance(self._device_num, Signature):
            device_num = self._device_num.var_name
        else:
            device_num = self._device_num

        # ok return it
        return f"acc set device_num({device_num})"

class ACCSetDeviceNumDirective(ACCDirective, RegionDirective, metaclass=abc.ABCMeta):
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
    def __init__(self, parent=None, children=None, device_num=None):
        # A Directive always contains a Schedule
        super().__init__(children=children, parent=parent)
        self._device_num = device_num

    def __eq__(self, other):
        '''
        Checks whether two nodes are equal. Two ACCKernelsDirective nodes are
        equal if their default_present members are equal.

        :param object other: the object to check equality to.

        :returns: whether other is equal to self.
        :rtype: bool
        '''
        is_eq = super().__eq__(other)
        is_eq = is_eq and self._device_num == other._device_num

        return is_eq

    @property
    def device_num(self):
        '''
        :returns: whether the "default(present)" clause is added to the \
                  kernels directive.
        :rtype: bool
        '''
        return self._device_num

    def gen_code(self, parent):
        '''
        Generate the f2pygen AST entries in the Schedule for this
        OpenACC Kernels directive.

        :param parent: the parent Node in the Schedule to which to add this \
                       content.
        :type parent: sub-class of :py:class:`psyclone.f2pygen.BaseGen`

        '''
        self.validate_global_constraints()

        # We re-use the 'begin_string' method but must skip the leading 'acc'
        # that it includes.
        parent.add(DirectiveGen(parent, "acc", "begin",
                                *self.begin_string().split()[1:]))

        self.gen_post_region_code(parent)

    def begin_string(self):
        '''Returns the beginning statement of this directive, i.e.
        "acc kernels ...". The backend is responsible for adding the
        correct directive beginning (e.g. "!$").

        :returns: the beginning statement for this directive.
        :rtype: str

        '''
        result = f"acc set device_num({self._device_num})"

        return result

    def end_string(self):
        '''
        Returns the ending statement for this directive. The backend is
        responsible for adding the language-specific syntax that marks this
        as a directive.

        :returns: the closing statement for this directive.
        :rtype: str

        '''
        return ""

class CrocoACCRaiseKernelThroughIf(Transformation):
    def __str__(self):
        return "Makes ACC kernel directive raising up through if/else"

    @property
    def name(self):
        return "CrocoACCRaiseKernelThroughIf"

    def apply(self, routine: Routine, options=None):
        # Ensure that the proposed transformation is valid
        self.validate(routine, options)

        # search all if/else
        if_blocks = routine.walk(IfBlock)
        for if_block in if_blocks:
            childs = 0
            leaf_child_is_uniq_acc_kernel = 0
            for branch in if_block.children:
                if isinstance(branch, Schedule):
                    childs += 1
                    if len(branch.children) == 1 and isinstance(branch.children[0], ACCKernelsDirective):
                        leaf_child_is_uniq_acc_kernel += 1
                    elif len(branch.children) == 0:
                        leaf_child_is_uniq_acc_kernel += 1

            print(if_block.view())
            print([leaf_child_is_uniq_acc_kernel, childs])
            if leaf_child_is_uniq_acc_kernel == childs and leaf_child_is_uniq_acc_kernel != 0:
                # Create a directive containing the nodes in node_list and insert it.
                parent = if_block.parent
                start_index = if_block.position
                directive = PsycloneACCKernelsDirective(
                    parent=parent, children=[if_block.detach()],
                    default_present=True, default_async='1')
                parent.children.insert(start_index, directive)

                # delete kernel directive in childs
                print(if_block.view())
                for branch in if_block.children:
                    # move up acc childs
                    if isinstance(branch, Schedule) and len(branch.children) > 0:
                        acc = branch.children[0]
                        assert isinstance(acc, PsycloneACCKernelsDirective)
                        assert isinstance(acc.children[0], Schedule)
                        for node in acc.children[0].children:
                            node.detach()
                            branch.children.append(node)

                        # remove acc kernek in if block
                        acc.detach()

    def validate(self, sched, options=None):
        pass

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

from poseidon.dsl.helper import *
from poseidon.dsl.kernel import PsycloneACCKernelsDirective
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.ifblock import IfBlock
from psyclone.psyir.nodes.call import Call
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.nodes import ACCDirective, RegionDirective
from psyclone.f2pygen import DirectiveGen
from poseidon_extensions import ACCWaitDirective, ACCSetDeviceNumDirective
from psyclone.transformations import ACCEnterDataTrans
from psyclone.psyGen import Transformation
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.nodes import ACCDataDirective, ACCDirective, \
    ACCEnterDataDirective, ACCKernelsDirective, ACCLoopDirective, \
    ACCParallelDirective, ACCRoutineDirective, Assignment, CodeBlock, \
    Directive, Loop, Node, OMPDeclareTargetDirective, \
    OMPDirective, OMPMasterDirective, \
    OMPParallelDirective, OMPParallelDoDirective, OMPSerialDirective, \
    OMPSingleDirective, OMPTaskloopDirective, PSyDataNode, Reference, \
    Return, Routine, Schedule, ACCUpdateDirective
from utils import add_kernels, normalise_loops, \
    insert_explicit_loop_parallelism
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from psyclone.psyir.transformations import ACCUpdateTrans
from poseidon_extensions import CrocoACCRaiseKernelThroughIf

from psyclone.core import Signature
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir

def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    for invoke in psy.invokes.invoke_list:
        # Convert array and range syntax to explicit loops
        normalise_loops(
            invoke.schedule,
            unwrap_array_ranges=True,
            hoist_expressions=False,
        )

    kernels = extract_kernels_from_psyir(psy.container)
    kernels.make_acc_tranformation(False)
    #kernels.merge_joinable_kernels()

    for kernel in kernels.kernels:
        loop = kernel.root_node

        
        try:
            ACCLoopTrans().apply(loop)
        except TransformationError as err:
            # This loop can not be transformed, proceed to next loop
            print("Loop not parallelised because:", str(err))
            continue

        # Count the number of perfectly nested loops & collapse them
        num_nested_loops = 0
        next_loop = loop
        while isinstance(next_loop, Loop):
            num_nested_loops += 1
            if len(next_loop.loop_body.children) > 1:
                break
            next_loop = next_loop.loop_body.children[0]

        if num_nested_loops > 1:
            loop.parent.parent.collapse = num_nested_loops

    routines = psy.container.walk(Routine)
    for routine in routines:
        if routine.name == "step2D_FB_tile":
            ifs = psy.container.walk(IfBlock)
            writer = FortranWriter ()
            for last_if in reversed(ifs):
                if writer(last_if.condition) == "iif == nfast":
                    last_if.if_body.children.insert(0, ACCWaitDirective(wait_queue=1))
                    break
        elif routine.name == "step2d":
            calls = psy.container.walk(Call)
            call = calls[0]
            pos = call.position
            call.parent.children.insert(pos, ACCSetDeviceNumDirective(device_num='tile'))

        #################################################################
        vars_to_fetch = []
        for ref in invoke.schedule.walk(Reference):
            if ref.name.lower() in ['dc1d','cd1d','ffc1d','fc1d','cf1d','dr1d','dz1d','bc1d']:
                vars_to_fetch.append(ref.name.lower())
        vars_to_fetch = list(set(vars_to_fetch))

        if len(vars_to_fetch) > 100000:
            # build sigs
            signatures = []
            for var in vars_to_fetch:
                print(var)
                signatures.append(Signature(var))
            print(signatures)

            ACCEnterDataTrans().apply(invoke.schedule, options={'signatures': signatures})
            for dir in invoke.schedule.walk(AccEnterDataDir):
                dir.sig_set = signatures

        #################################################################
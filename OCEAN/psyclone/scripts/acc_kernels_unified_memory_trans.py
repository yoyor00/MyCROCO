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

from utils import add_kernels

from psyclone.psyir.nodes import Directive, Reference, ACCUpdateDirective
from psyclone.psyir.transformations import ACCUpdateTrans
from psyclone.transformations import ACCEnterDataTrans
from psyclone.core import Signature

from psyclone.core import Signature
from psyclone.nemo import NemoACCEnterDataDirective as \
                AccEnterDataDir


def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    print("Invokes found:")
    print("\n".join([str(name) for name in psy.invokes.names]))

    for invoke in psy.invokes.invoke_list:

        if not invoke.schedule:
            print(f"Invoke {invoke.name} has no Schedule! Skipping...")
            continue

        add_kernels(invoke.schedule.children, default_present=False)

        #################################################################
        vars_to_fetch = []
        for ref in invoke.schedule.walk(Reference):
            if ref.name.lower() in ['dc1d','cd1d','ffc1d','fc1d','cf1d','dr1d','dz1d','bc1d']:
                vars_to_fetch.append(ref.name.lower())
        vars_to_fetch = list(set(vars_to_fetch))

        if len(vars_to_fetch) > 0:
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

        
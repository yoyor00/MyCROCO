##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

'''A transformation script that seeks to apply OpenACC KERNELS directives to
NEMO style code. In order to use it you must first install PSyclone. See
README.md in the top-level directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s <this_script> <target_source_file>

The transformation script attempts to insert Kernels directives at the
highest possible location(s) in the schedule tree (i.e. to enclose as
much code as possible in each Kernels region).

'''


##########################################################
# python
import os
# poesidon
from lib.poseidon.dsl.helper import extract_kernels_from_psyir
from lib.extensions import loops
from lib.croco import psyacc
# intenral
from lib.extensions.directives import ACCWaitDirective, ACCSetDeviceNumDirective
from lib.from_nemo_psyclone.utils import normalise_loops
# psyclone
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes.if_block import IfBlock
from psyclone.psyir.nodes.call import Call
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.psyir.transformations.transformation_error import \
    TransformationError
from psyclone.psyir.nodes import Loop, Routine, IntrinsicCall
from psyclone.transformations import ACCLoopTrans

##########################################################
def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    # enable dump
    dump_snippets = False
    if 'CROCO_PSYCLONE_SNIPPET_DUMPS' in os.environ and os.environ['CROCO_PSYCLONE_SNIPPET_DUMPS'] == 'true':
        loops.ENABLE_SNIPPET_DUMPS = True
        dump_snippets = True

    # get routine walker from 'psy' ir
    routines = psy.container.walk(Routine)
    for routine in routines:
        psyacc.apply_psyclone_original_trans(routine, dump_snippets = dump_snippets)

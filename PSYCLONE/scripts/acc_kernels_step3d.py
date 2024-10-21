##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

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
# specific CROCO vars
VARS_1D = ['fc', 'cf', 'dc', 'bc', 'dz', 'dr']
VARS_3D = ['fx', 'fe', 'work', 'work2']

##########################################################
# python
import os
# internal
from lib.extensions import acc
from lib.extensions import kernels
from lib.extensions import scratch
from lib.extensions import loops
from lib.croco import step3d
# poseidon
from lib.poseidon.dsl.helper import *
from lib.extensions import loops
# psyclone
from psyclone.psyir.nodes.routine import Routine
from psyclone.psyir.nodes import Loop, Routine

##########################################################
def trans(psy):
    """
    A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' to Croco code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    """

    # enable dump
    dump_snippets = False
    if 'CROCO_PSYCLONE_SNIPPET_DUMPS' in os.environ and os.environ['CROCO_PSYCLONE_SNIPPET_DUMPS'] == 'true':
        dump_snippets = True

    # steps
    acc.add_missing_device_vars(psy.container)
    acc.set_device_tile(psy.container)

    routines = psy.container.walk(Routine)
    for routine in routines:  
        step3d.apply_step3d_routine_trans(routine, psy.container, dump_snippets = dump_snippets)

    ##################################################################
    # automatic bench to fix the private issue
    #apply_acc_fetch_vars(psy)

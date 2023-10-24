##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
#PSyIR internal representation
from psyclone.psyir.nodes import Loop, ArrayReference
# internal
from . walker import WalkerCallbackInterface
from ..dsl.kernel import Kernel, KernelList
from ..base.types import AccessMode

##########################################################
class KernelExtractor(WalkerCallbackInterface):
    def __init__(self) -> None:
        self.loop_depth = 0
        self.kernel = None
        self.kernels = KernelList()

    def on_loop(self, loop: Loop, enter: bool) -> None:
        if enter:
            if self.kernel == None:
                assert self.loop_depth == 0
                self.kernel = Kernel(loop)
            self.kernel.push_loop_var(loop.variable.name)
            self.loop_depth += 1
        else:
            self.loop_depth -= 1
            if self.loop_depth == 0 and self.kernel != None:
                self.kernels.push_kernel(self.kernel)
                self.kernel = None

    def on_array_reference(self, array: ArrayReference, access_mode: AccessMode, enter: bool) -> None:
        if self.kernel != None:
            self.kernel.analyse_stencil_mask(array, access_mode)
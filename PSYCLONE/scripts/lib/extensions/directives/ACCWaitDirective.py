##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# psyclone
from psyclone.core import Signature
from psyclone.psyir.nodes import ACCStandaloneDirective
from psyclone.f2pygen import DirectiveGen

##########################################################
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

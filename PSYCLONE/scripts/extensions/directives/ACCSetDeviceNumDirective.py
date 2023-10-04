##########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement some new type of ACC nodes to help porting Croco.
(work ongoing to merge them in official psyclone).
'''

##########################################################
from psyclone.psyir.nodes import ACCDirective, RegionDirective
from psyclone.f2pygen import DirectiveGen
import abc
from psyclone.psyir.nodes import ACCDirective

##########################################################
class ACCSetDeviceNumDirective(ACCDirective, RegionDirective, metaclass=abc.ABCMeta):
    '''
    Represent the "acc set devicenum()" directive.
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

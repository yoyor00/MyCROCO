##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
from abc import ABC, abstractmethod
# internal
from ..config import Config

##########################################################
class AbstractCrocoSetup(ABC):
    '''
    Base class to be extended to support the new MiniCroco with CMake of the
    old default build system to run on official master branch.
    '''

    def __init__(self, config: Config, builddir: str):
        '''
        Constructor of the croco set.

        Parameters
        ----------
        config: Config
            To get access to the global config if needed.
        builddir: str
            Path of the build directory we want to create and build in.
        '''
        self.config = config
        self.sourcedir = config.croco_source_dir
        self.builddir = os.path.abspath(builddir)

    @abstractmethod
    def configure(self, configure_args: str) -> None:
        '''
        Perform the configuration of CROCO to be ready to build.

        Parameters
        ----------
        configure_args: str
            The arguments as a string as they would have been given to the configure script wrapping cmake.
            For the old corco variant it will automatically be parsed and transformed for compatibility
            translation.
        '''
        raise Exception("Not provided by implementation !")

    @abstractmethod
    def make(self, make_jobs: str) -> None:
        '''
        Perform the compilation either with jobcomp of make.

        Parameters
        ----------
        make_jobs: str
            To pass the `-j8` option if wanted.
        '''
        raise Exception("Not provided by implementation !")

    @abstractmethod
    def convert_patch_fnames(self, filename: str) -> str:
        '''
        Some files have been renamed between minicroco and the original one,
        in order to keep the whole case definition in bench agnostic to this
        we rename on the fly some files.
        '''
        raise Exception("Not provided by implementation !")

    @abstractmethod
    def copy_config(self, refdir_case: str, case_name: str, case_patches: dict) -> None:
        '''
        Copy the required config files to pass the case in the reference dir so
        we can possibly also reproduce by hand easily if needed one day
        (outside of BENCH).

        Parameters
        ----------
        refdir_case: str
            The path in which to put the files.
        case_name: str
            Name of the case to know which file to take (possibly).
        case_patches: dict
            The information from the case on what was patched.
        '''
        raise Exception("Not provided by implementation !")

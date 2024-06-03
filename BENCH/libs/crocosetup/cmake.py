##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json
import shutil
# internal
from ..helpers import Messaging, patch_lines, move_in_dir, run_shell_command
from .setup import AbstractCrocoSetup
from ..config import Config

##########################################################
class CMakeCrocoSetup(AbstractCrocoSetup):
    '''
    Base class to be extended to support the new MiniCroco with CMake of the
    old default build system to run on official master branch.
    '''

    def __init__(self, config: Config, builddir: str):
        super().__init__(config, builddir)

    def configure(self, configure_command: str) -> None:
        with move_in_dir(self.builddir):
            run_shell_command(configure_command, capture=self.config.capture)

    def cppdef_h_set_key(self, key_name: str, status: bool):
        '''
        Force enabling or disabling some cppkeys in cppdefs.h. In order to
        be applied to all cases without having to search the in case position.

        We place it just after the cases and before inclusion of cppdefs_dev.h
        which need to have everything well defined.

        Note: It is necessary to put # undef before # define as some keys can be
              already defined.
        '''
        # caution here, the space in '# define' is required due to croco tricks
        # do not remove it or it behave strangely
        if status:
            to_insert = f'# undef {key_name}\n# define {key_name}\n'
        else:
            to_insert = f'# undef {key_name}\n'

        # build the patch rule
        rules = [{
            'mode'  : 'insert-at-end',
            'what': to_insert,
            'descr' : f"Set {key_name} to {status}"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, 'cppdefs_override.h'), rules)

    def make(self, make_jobs: str) -> None:
        with move_in_dir(self.builddir):
            run_shell_command(f"make {make_jobs}", capture=self.config.capture)

    def convert_patch_fnames(self, filename: str) -> str:
        '''
        Some files have been renamed between minicroco and the original one,
        in order to keep the whole case definition in bench agnostic to this
        we rename on the fly some files.
        '''
        if filename == 'param.h':
            return 'param_override.h'
        else:
            return filename

    def copy_config(self, refdir_case: str, case_name: str, patches_and_keys: dict) -> None:
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
        patches_and_keys: dict
            The information from the case on what was patched.
        '''

        # vars
        builddir = self.builddir

        # create subdir
        put_int = f"{refdir_case}/cmake-mode/"
        os.makedirs(put_int, exist_ok=True)
        os.makedirs(f"{put_int}/OCEAN", exist_ok=True)

        # copy files from build sys
        to_copy = [
            # to know which cmake command has been used
            "configure.log",
            "cppdefs_override.h",
            "cppdefs_dev_override.h",
            "param_override.h",
            "OCEAN/config.h",
            "OCEAN/config_post.h",
        ]

        # copy them
        for file in to_copy:
            shutil.copyfile(f"{builddir}/{file}", f"{put_int}/{file}")

        # write extra
        with open(f"{put_int}/patches-and-keys.json", "w+") as fp:
            json.dump({
                "case": case_name,
                "changes": patches_and_keys
            }, fp)

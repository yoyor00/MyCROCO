##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import shlex
import argparse
import shutil
# internal
from ..helpers import Messaging, patch_lines, move_in_dir, run_shell_command
from .setup import AbstractCrocoSetup
from ..config import Config

##########################################################
class JobcompCrocoConfig:
    '''
    In the old croco there is several tricks to setup OpenMP.....
    We agregate here them to keep OldCrocoSetup readable.
    '''

    def __init__(self, builddir: str, minicroco: bool = False):
        self.builddir = builddir
        self.minicroco = minicroco

    def cppdef_h_select_case(self, case_name: str):
        '''
        Patch the cppdef.h to select the CASE to be used.
        '''

        # progress
        Messaging.step(f"Select case : {case_name}")

        # set patching rules
        rules = []

        # disable REGINAL (enabled by default in CROCO master, disabled by default in minicroco branch)
        if not self.minicroco:
            rules.append({
                'mode'  : 'replace', 
                'what'  : "#define REGIONAL        /* REGIONAL Applications */\n", 
                'by'    : "#undef REGIONAL        /* REGIONAL Applications */\n",
                'descr' : 'Disable default REGIONAL case'
            })

        # insert the wated one
        rules.append({
            'mode'  : 'insert-after', 
            'what'  : "#undef REGIONAL        /* REGIONAL Applications */\n", 
            'insert': f'#define {case_name}\n',
            'descr' : f'Enabled wanted {case_name} case'
        })

        # apply
        patch_lines(os.path.join(self.builddir, 'cppdefs.h'), rules)

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
            'mode'  : 'insert-after-before',
            'after' : f'#endif /* END OF CONFIGURATION CHOICE */', 
            'before': '#include "cppdefs_dev.h"',
            'insert': to_insert,
            'descr' : f"Set {key_name} to {status}"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, 'cppdefs.h'), rules)

    def param_h_configure_openmp_split(self, threads:int):
        '''
        Configure the OpenMP splitting in param.h
        '''

        # build the patching rule
        rules = [{
            'mode'  : 'replace', 
            'what'  :  "      parameter (NPP=4)\n", 
            'by'    : f"      parameter (NPP={threads})\n",
            'descr' : f"Set OPENMP threads = {threads}"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, 'param.h'), rules)

    def jobcomp_configure_set_extra_fflags(self, fflags: str):
        '''
        Permit to append some extra fortran flags for the compiler. Typically
        the `-march` one to be tuned for the local CPU we want.
        '''

        # build the rule
        rules = [{
            'mode'  : 'insert-after-before',
            'after' : 'FFLAGS1=${CROCO_FFLAGS1-$FFLAGS1}',
            'before': 'CFT1=${CROCO_CFT1-$CFT1}',
            'insert': f'FFLAGS1=$FFLAGS1 {fflags}\n',
            'descr' : f"Set FFLAGS to {fflags}"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, 'jobcomp'), rules)

    def jobcomp_configure_set_compiler(self, fortran_compiler: str):
        '''
        Change the compiler.
        '''

        # build the rull
        if self.minicroco:
            rules = [{
                'mode'  :  'replace',
                'what'  :  'test -z "$!FC" && FC=gfortran\n',
                'by'    : f'test -z "$!FC" && FC={fortran_compiler}\n',
                'descr' : f"Set fortran compiler to FC={fortran_compiler}"
            }]
        else:
            rules = [{
                'mode'  :  'replace',
                'what'  :  'FC=gfortran\n',
                'by'    : f'FC={fortran_compiler}\n',
                'descr' : f"Set fortran compiler to FC={fortran_compiler}"
            }]

        # apply
        patch_lines(os.path.join(self.builddir, 'jobcomp'), rules)

    def param_h_configure_mpi_split(self, splitting: str):
        '''
        Configure the MPI domain splitting and usable number of ranks.

        Parameters
        ----------
        splitting: str
            Splitting under the form of 2x8 as used for the CMake scripts.
        '''
        # split
        dims = splitting.split('x')
        assert len(dims) == 2

        # extrat
        np_x = int(dims[0])
        np_eta = int(dims[1])
        ranks = np_x * np_eta

        # build rule
        rules = [{
            'mode'  : 'replace', 
            'what'  :  "      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)\n", 
            'by'    : f"      parameter (NP_XI={np_x},  NP_ETA={np_eta},  NNODES=NP_XI*NP_ETA)\n",
            'descr' : f"Set MPI splitting : {splitting} : ranks={ranks}, np_xi={np_x}, np_eta={np_eta}"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, 'param.h'), rules)

    def fix_jobcomp_missing_return_status(self) -> None:
        '''
        Some fixes (to be merged once).

        This one fix the fact that if `make` fails, jobcomp currently still
        return an OK status.

        TODO: probably also make sure that the `croco` binary is not there anymore
        if build has failed. Not using it without noticing it has not been rebuilt.
        '''

        # build rule
        rules = [{
            'mode'  : 'replace', 
            'what'  :  "$MAKE\n", 
            'by'    : f"$MAKE || exit 1\n",
            'descr' : f"Fix missing return make status in jobcomp"
        }]

        # apply
        patch_lines(os.path.join(self.builddir, "jobcomp"), rules, allow_already_done=True)

##########################################################
class JobcompCrocoSetup(AbstractCrocoSetup):
    '''
    Implement configure wrapping of old croco build system to be able to also
    use bench on the master branch of croco.
    '''

    def __init__(self, config: Config, builddir: str):
        super().__init__(config, builddir)
        self.croco_config = JobcompCrocoConfig(builddir, minicroco=config.is_minicroco_jobcomp)

    @staticmethod
    def convert_arg_for_argparse(arg_string: str) -> list:
        '''
        Splits the arguments as passed to cmake to get them as a list of args to
        be passed to argparser.

        Note: We need to reshape a bit, because the configure script coming with
              cmake is using `--with-xxx=value`, where argparse want with space :
              `--with-xxx value`.
        '''

        # make splitting by keeping the spaces protected by "" or '' as in shell.
        basic_list = shlex.split(arg_string)

        # need to make vars (FC=gfortran) close together because argparse don't
        # like to have them mixed with options
        tmp_opts = []
        tmp_vars = []

        # pass over the args
        for entry in basic_list:
            # need to split --with-xxx=value => ['--with-xxx', 'value']
            if entry.startswith("--") and '=' in entry:
                tmp_opts += entry.split('=', maxsplit=1)
            else:
                tmp_vars.append(entry)

        # fuse again vars + args
        return tmp_vars[1:] + tmp_opts

    def cppdef_h_set_openmp(self, status: bool):
        '''
        Enable of disable OpenMP
        '''
        self.croco_config.cppdef_h_set_key('OPENMP', status)

    def cppdef_h_enable_openacc(self, status: bool):
        '''
        Enable OPENACC
        '''
        self.croco_config.cppdef_h_set_key('OPENACC', status)

    def cppdef_h_enable_openacc_psyclone(self, status: bool):
        '''
        Enable OPENACC_PSYCLONE
        '''
        self.croco_config.cppdef_h_set_key('OPENACC_PSYCLONE', status)

    def convert_patch_fnames(self, filename: str) -> str:
        '''
        Some files have been renamed between minicroco and the original one,
        in order to keep the whole case definition in bench agnostic to this
        we rename on the fly some files.
        '''
        if filename == 'param_override.h':
            return 'param.h'
        else:
            return filename

    def handle_variables(self, arg_vars: list, is_not_mpi: bool) -> dict:
        '''
        Apply to ops required when getting some variables on the command line.
        
        Currently consider (exemple of values): 
         - FC=gfortran
         - FFLAGS=-march=native
        '''
        vars = {}
        for entry in arg_vars:
            # split name & value
            fields = entry.split('=', maxsplit=1)

            # extract
            var_name = fields[0]
            var_value = fields[1]

            # store if we want to use latter
            vars[var_name] = var_value

            # apply vars
            if var_name == 'FFLAGS':
                self.croco_config.jobcomp_configure_set_extra_fflags(var_value)
            elif var_name == 'FC':
                # there is an issue if we force mpifort that way with jobcomp
                if is_not_mpi:
                    self.croco_config.jobcomp_configure_set_compiler(var_value)
            else:
                raise Exception(f"Unsupported variable : {entry}")

        #ok
        return vars

    def emulate_cmd_line_parse_args(self, minicroco_args: str):
        '''
        Convert the original CMake + configure semantic in the old one.

        Note: This might permit to one day use this script as a configure to
              be added in the old CROCO.
        '''

        # convert conventional configure agument way (with --xxx=v => -xxx v)
        args = self.convert_arg_for_argparse(minicroco_args)

        # build parser
        parser = argparse.ArgumentParser(
            prog = 'old_croco_configure',
            description = 'A wrapper to amulate the new build API onto the old croco.')
        
        # add options
        parser.add_argument("--with-case", required=True, help="Select the case to run.")
        parser.add_argument("--with-optim", required=False, help="Select the optimization mode to use.", default='seq')
        parser.add_argument("--with-threads", required=False, help="Select the number of threads.", default='4')
        parser.add_argument("--with-splitting", required=False, help="Select the MPI domain splitting.", default='1x4')
        parser.add_argument("VARS", nargs='*', type=str, help="Extra variable definitions, like compilers : FC=gfortran.")

        # parser
        options = parser.parse_args(args)

        # ok
        return options

    def create_build_dir(self) -> None:
        '''
        Create the build directory by calling `create_config.bash` so it copies
        everything.
        '''
        # create directory
        with move_in_dir(self.sourcedir):
            # because this !**ù$$m script does not like getting a full path
            rel_path = os.path.relpath(self.builddir)
            run_shell_command(f"./create_config.bash -f -n {rel_path}", capture=self.config.capture)

    def configure(self, minicroco_args: str) -> None:
        '''
        Perform the configuration.

        Parameters
        ----------
        minicroco_args: str
            Arguments on the form noramly passed to the minicroco CMake + configure
            script which will be translated in what is needed for the old CROCO.
        '''

        # configure
        Messaging.step(f"jobcomp-configure {minicroco_args}")

        # get options
        options = self.emulate_cmd_line_parse_args(minicroco_args)

        # create build dir
        self.create_build_dir()

        # set case
        self.croco_config.cppdef_h_select_case(options.with_case)

        # config
        use_openmp = False
        use_mpi = False
        use_openacc = False
        use_openacc_psyclone = False

        # apply variable
        is_not_mpi = (options.with_optim != 'mpi')
        vars = self.handle_variables(options.VARS, is_not_mpi)

        # apply optim mode
        if options.with_optim == 'seq':
            pass
        elif options.with_optim == 'openmp':
            use_openmp = True
        elif options.with_optim == 'mpi':
            use_mpi = True
        elif options.with_optim == 'openacc-native':
            use_openacc = True
            use_openacc_psyclone = False
        elif options.with_optim == 'openacc-psyclone':
            use_openacc = True
            use_openacc_psyclone = True
        elif options.with_optim == 'poseidon':
            raise Exception(f"Variant {options.with_optim} not yet supported")
        else:
            raise Exception("Invalid optimisation mode from commande line : {minicroco_args}")
        
        # set config
        self.cppdef_h_set_openmp(use_openmp)
        self.cppdef_h_enable_openacc(use_openacc)
        self.cppdef_h_enable_openacc_psyclone(use_openacc_psyclone)
        if use_openmp:
            self.croco_config.param_h_configure_openmp_split(options.with_threads)
        if use_mpi:
            self.croco_config.param_h_configure_mpi_split(options.with_splitting)

        # fix issue
        self.croco_config.fix_jobcomp_missing_return_status()

    def make(self, make_jobs: str):
        '''
        Perform the build with jobcomp.
        '''

        # to pass -j8 we need to go through env var
        os.environ["MAKEFLAGS"] = make_jobs

        # move in dir & call jobcomp
        with move_in_dir(self.builddir):
            run_shell_command(f"./jobcomp", capture=self.config.capture)

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

        # vars
        builddir = self.builddir

        # create subdir
        put_int = f"{refdir_case}/jobcomp-mode/"
        os.makedirs(put_int, exist_ok=True)

        # copy files from build sys
        to_copy = [
            # keep track of the exact case config being used.
            "cppdefs.h",
            "cppdefs_dev.h",
            "param.h",
        ]

        # copy them
        for file in to_copy:
            shutil.copyfile(f"{builddir}/{file}", f"{put_int}/{file}")

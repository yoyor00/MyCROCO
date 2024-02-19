##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json
import shutil
import platform
import subprocess
from .config import Config
from .helpers import move_in_dir, run_shell_command, apply_vars_in_str, Messaging, patch_lines
from .hyperfine import run_hyperfine
from .check import compare_netcdf_files
from .crocosetup.cmake import CMakeCrocoSetup
from .crocosetup.jobcomp import JobcompCrocoSetup

##########################################################
class Croco:
    def __init__(self, config: Config, case_name: str, variant_name: str):
        # keep track
        self.config = config
        self.case_name = case_name
        self.variant_name = variant_name
        self.full_name = f"{case_name}-{variant_name}"

        # build directory name
        self.dirname = f"{config.workdir}/{case_name}/{variant_name}"
        self.dirname = os.path.abspath(self.dirname)

        # extract infos
        self.case = config.config['cases'][case_name]
        self.variant = config.config['variants'][variant_name]

        # if use old croco
        #self.old_croco = None
        if config.has_cmake:
            self.croco_build = CMakeCrocoSetup(self.config, self.dirname)
        else:
            self.croco_build = JobcompCrocoSetup(self.config, self.dirname)

    def reset(self):
        # display
        Messaging.step("Reset")

        # rm the old one
        shutil.rmtree(self.dirname, ignore_errors=True)

    def configure(self):
        # Configure
        Messaging.step("Configure")

        # extract some needed vars
        croco_source_dir = self.config.croco_source_dir
        dirname = self.dirname
        case_cpp_key = self.case['case']
        tuning_familly = self.variant['tuning_familly']
        tuning_flags = self.config.host['tuning'][tuning_familly]
        croco_build = self.croco_build

        # create dir
        os.makedirs(dirname, exist_ok=True)

        # extract options
        configure_variant_options = self.variant['configure']
        configure_case_option = f"--with-case={case_cpp_key}"
        if tuning_flags != '':
            configure_compiler_option = f"FFLAGS=\"{tuning_flags}\""
        else:
            configure_compiler_option = ""

        # build command
        command = f"{croco_source_dir}/configure {configure_compiler_option} {configure_case_option} {configure_variant_options}"
        command = apply_vars_in_str(command, {'case': self.case, 'tuning': self.config.host['tuning']})

        # jump in & configure
        croco_build.configure(command)

    def compile(self):
        # display
        Messaging.step("Compile")

        # extract some needed vars
        croco_build = self.croco_build
        make_jobs = f"-j{self.config.make_jobs}"

        # jump in & build
        croco_build.make(make_jobs)

    def build(self, extra_info: str = "", force_rebuild: bool = False):
        # display
        Messaging.section(f"Building CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")
        # perform steps
        if force_rebuild:
            self.reset()
        elif not force_rebuild and os.path.exists(f"{self.dirname}/croco"):
            Messaging.step(f"Already built...")
            return

        # effectively built
        self.configure()
        self.setup_case()
        self.compile()

    def run(self, extra_info: str = ""):
        # display
        Messaging.section(f"Running CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")

        # extract some vars
        command_prefix = self.variant['command_prefix']
        environ = self.variant.get('environ', {}).copy()
        dirname = self.dirname
        runs = self.config.runs
        results = self.config.results
        host_tuning = self.config.host['tuning']

        # load tuning env & override if needed
        for key, value in host_tuning['environ'].items():
            environ[key] = value

        # apply vars
        command_prefix = apply_vars_in_str(command_prefix, {'case': self.case, 'tuning': host_tuning})

        # build end
        env_line = ""
        for var, value in environ.items():
            env_line += f"{var}=\"{value}\" "

        # build command and run
        command = f"{env_line} ../../../scripts/correct_end.sh {command_prefix} ./croco"
        with move_in_dir(dirname):
            result = run_hyperfine(command, runs=runs, verbose=self.config.verbose)

        # calc output
        output_filename = f"{results}/result-{self.full_name}.json"

        # add some infos
        result['name'] = self.full_name
        result['original_path'] = output_filename
        result['hostname'] = platform.node()

        # dump results
        with open(output_filename, "w+") as fp:
            json.dump(result, fp=fp, indent='\t')

    def setup_case(self):
        # apply the case paches
        case_name = self.case_name
        patches = self.case['patches']
        croco_build = self.croco_build

        # display
        Messaging.step(f"Apply case config : {case_name}")

        # loop
        with move_in_dir(self.dirname):
            for file, changes in patches.items():
                # convert to list if needed
                if isinstance(changes, dict):
                    changes = [changes]

                # change file name for old croco
                file_filtered = croco_build.convert_patch_fnames(file)

                # loop on all changes
                for change in changes:
                    patch_lines(file_filtered, [
                        {
                            "mode": "replace",
                            "after": change['next_to'],
                            "what": change['what'] + '\n',
                            "by": change['by'] + '\n',
                            "descr": change['descr']
                        }
                    ])

    def check_one_file_from_seq_ref(self, filename: str) -> None:
        # extract vars
        case_name = self.case_name
        dirname = self.dirname

        # if sequential skip
        if self.variant_name == 'sequential':
            Messaging.step(f"Checking {case_name} / {filename} skiped for sequential")
            return
        else:
            Messaging.step(f"Checking {case_name} / {filename}")

        # error
        seq_file = f"{dirname}/../sequential/{filename}"
        if not os.path.exists(seq_file):
            raise Exception(f"Missing '{seq_file}', are you sure you ran case 'sequential' first to get a reference for checks ?")

        # compare
        actual_file = f"{dirname}/{filename}"
        compare_netcdf_files(seq_file, actual_file)

    def check_one_file_from_ref_dir(self, filename: str) -> None:
        # vars
        refdir = self.config.use_ref
        case_name = self.case_name
        dirname = self.dirname

        # progress
        Messaging.step(f"Checking {case_name} / {filename} from ref={refdir}")

        # error
        ref_file = f"{refdir}/{case_name}/{filename}"
        if not os.path.exists(ref_file):
            raise Exception(f"Missing '{ref_file}', are you sure you ran case {case_name} in reference directory {refdir} ?")
        
        # compare
        actual_file = f"{dirname}/{filename}"
        compare_netcdf_files(ref_file, actual_file)

    def check_one_file(self, filename: str):
        # vars
        refdir = self.config.use_ref

        # swtich
        if refdir:
            self.check_one_file_from_ref_dir(filename)
        else:
            self.check_one_file_from_seq_ref(filename)

    def check(self):
        for filename in self.case['check_outputs']:
            self.check_one_file(filename)

    def make_ref(self):
        # extract vars
        refdir = self.config.build_ref
        case_name = self.case_name
        dirname = self.dirname
        case_patches = self.case['patches']

        # check
        assert refdir

        # progress
        Messaging.step(f"Storing ref {case_name} into ref={refdir}")

        # we create only for seq
        if self.variant_name != 'sequential':
            return

        # create dir if not exist
        refdir_case = f"{refdir}/{case_name}"
        os.makedirs(refdir_case, exist_ok=True)

        # copy the files there
        for filename in self.case['check_outputs']:
            orig_path = f"{dirname}/{filename}"
            dest_path = f"{refdir_case}/{filename}"
            Messaging.step(f"Storing file {case_name} / {filename}")
            shutil.copyfile(orig_path, dest_path)

        # also copy the config
        self.croco_build.copy_config(refdir_case, case_name, case_patches)

        # add the case config files
        case_capitalized = case_name.capitalize()
        case_file = f"TEST_CASES/croco.in.{case_capitalized}"
        os.makedirs(f"{refdir}/{case_name}/TEST_CASES", exist_ok=True)
        shutil.copyfile(f"{dirname}/{case_file}", f"{refdir}/{case_name}/{case_file}")

        # copy the case file under croco.in
        if not os.path.exists(f"{refdir}/{case_name}/croco.in"):
            os.symlink(case_file, f"{refdir}/{case_name}/croco.in")

        # dump case info
        with open(f"{refdir}/{case_name}/case.json", "w+") as fp:
            json.dump(self.case, fp, indent='\t')

        # dump git info
        with open(f"{refdir}/{case_name}/git.log", "w+") as fp:
            fp.write(subprocess.getoutput("git log HEAD~1..HEAD 2>/dev/null || echo 'No git !'") + "\n")

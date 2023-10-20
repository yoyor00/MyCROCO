##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json
import shutil
from .config import Config
from .helpes import move_in_dir, run_shell_command, apply_vars_in_str, Messaging, patch_lines
from .hyperfine import run_hyperfine

##########################################################
class Croco:
    def __init__(self, config: Config, case_name: str, variant_name: str):
        # keep track
        self.config = config
        self.case_name = case_name
        self.variant_name = variant_name
        self.full_name = f"{case_name}-{variant_name}"

        # build directory name
        self.dirname = f"{config.workdir}/build-{self.full_name}"
        self.dirname = os.path.abspath(self.dirname)

        # extract infos
        self.case = config.config['cases'][case_name]
        self.variant = config.config['variants'][variant_name]

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
        verbose = self.config.verbose

        # create dir
        if not os.path.exists(dirname):
            os.mkdir(dirname)

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
        with move_in_dir(dirname):
            run_shell_command(command, capture=not verbose)

    def compile(self):
        # display
        Messaging.step("Compile")

        # extract some needed vars
        dirname = self.dirname
        make_jobs = f"-j{self.config.make_jobs}"

        # jump in & build
        with move_in_dir(dirname):
            run_shell_command(f"make {make_jobs}")

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
        self.compile()
        self.setup_case()

    def run(self, extra_info: str = ""):
        # display
        Messaging.section(f"Running CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")
        # extract some vars
        command_prefix = self.variant['command_prefix']
        environ = self.variant.get('environ', {})
        dirname = self.dirname
        runs = self.config.runs
        results = self.config.results

        # apply vars
        command_prefix = apply_vars_in_str(command_prefix, {'case': self.case, 'tuning': self.config.host['tuning']})

        # build end
        env_line = ""
        for var, value in environ.items():
            env_line += f"{var}=\"{value}\""

        # build command and run
        command = f"{env_line} ../../scripts/correct_end.sh {command_prefix} ./croco"
        with move_in_dir(dirname):
            result = run_hyperfine(command, runs=runs)

        # dump results
        with open(f"{results}/result-{self.full_name}.json", "w+") as fp:
            json.dump(result, fp=fp, indent='\t')

    def setup_case(self):
        # apply the case paches
        case_name = self.case_name
        patches = self.case['patches']

        # display
        Messaging.step(f"Apply case config : {case_name}")

        # loop
        with move_in_dir(self.dirname):
            for file, changes in patches.items():
                # convert to list if needed
                if isinstance(changes, dict):
                    changes = [changes]

                # loop on all changes
                for change in changes:
                    patch_lines(file, [
                        {
                            "mode": "replace",
                            "after": change['after'],
                            "what": change['what'] + '\n',
                            "by": change['by'] + '\n',
                        }
                    ])

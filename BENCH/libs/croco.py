##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import shutil
from .config import Config
from .helpes import move_in_dir, run_shell_command, apply_vars_in_str

##########################################################
class Croco:
    def __init__(self, config: Config, case_name: str, variant_name: str):
        # keep track
        self.config = config
        self.case_name = case_name
        self.variant_name = variant_name

        # build directory name
        self.dirname = f"{config.workdir}/build-{case_name}-{variant_name}"
        self.dirname = os.path.abspath(self.dirname)

        # extract infos
        self.case = config.config['cases'][case_name]
        self.variant = config.config['variants'][variant_name]

    def reset(self):
        # rm the old one
        shutil.rmtree(self.dirname, ignore_errors=True)

    def configure(self):
        # extract some needed vars
        croco_source_dir = self.config.croco_source_dir
        dirname = self.dirname
        case_cpp_key = self.case['case']

        # create dir
        if not os.path.exists(dirname):
            os.mkdir(dirname)

        # extract options
        configure_variant_options = self.variant['configure']
        configure_case_option = f"--with-case={case_cpp_key}"

        # build command
        command = f"{croco_source_dir}/configure {configure_case_option} {configure_variant_options}"
        command = apply_vars_in_str(command, {'case': self.case})

        # jump in & configure
        with move_in_dir(dirname):
            run_shell_command(command)

    def compile(self):
        # extract some needed vars
        dirname = self.dirname
        make_jobs = f"-j{self.config.make_jobs}"

        # jump in & build
        with move_in_dir(dirname):
            run_shell_command(f"make {make_jobs}")

    def run(self):
        self.reset()
        self.configure()
        self.compile()

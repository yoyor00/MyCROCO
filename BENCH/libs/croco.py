##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import glob
import json
import shutil
import platform
import subprocess

# internal
from .config import Config
from .helpers import (
    move_in_dir,
    run_shell_command,
    apply_vars_in_str,
    Messaging,
    patch_lines,
    copy_tree_with_absolute_symlinks,
)
from .hyperfine import run_hyperfine
from .check import compare_netcdf_files
from .crocosetup.cmake import CMakeCrocoSetup
from .crocosetup.jobcomp import JobcompCrocoSetup
from .meshdrawer import MeshDrawer


##########################################################
class Croco:
    def __init__(self, config: Config, case_name: str, variant_name: str):
        # keep track
        self.config = config
        self.case_name = case_name
        self.variant_name = variant_name
        self.full_name = f"{case_name}-{variant_name}"

        # build directory name
        self.dirname = self.calc_rundir(variant_name, case_name)
        self.dirname = os.path.abspath(self.dirname)
        self.dirname_result = self.calc_resdir(variant_name, case_name)
        self.dirname_result = os.path.abspath(self.dirname_result)

        # extract infos
        self.case = config.config["cases"][case_name]
        self.variant = config.config["variants"][variant_name]

        # if use old croco
        # self.old_croco = None
        if config.has_cmake and not self.config.force_jobcomp:
            self.croco_build = CMakeCrocoSetup(self.config, self.dirname)
        else:
            self.croco_build = JobcompCrocoSetup(
                self.config, self.dirname, self.variant["tuning_familly"]
            )

    def calc_rundir(self, variant_name, case_name):
        return f"{self.config.workdir}/{case_name}/{variant_name}"

    def calc_resdir(self, variant_name, case_name):
        return f"{self.config.results}/{case_name}/{variant_name}"

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
        dirname_result = self.dirname_result
        case_cpp_key = self.case["case"]
        tuning_familly = self.variant["tuning_familly"]
        tuning_fflags = self.variant["tuning_fflags"]
        tuning_flags = self.config.host["tuning"][tuning_familly][tuning_fflags]
        tuning_flags_extra = self.config.host["tuning"][tuning_familly]["extra"]
        croco_build = self.croco_build

        # create dir
        os.makedirs(dirname, exist_ok=True)
        os.makedirs(dirname_result, exist_ok=True)
        # link data dir
        if "croco_files_path" in self.case:
            self.input_dir = os.path.join(
                self.config.data_root_path, self.case["croco_files_path"]
            )
            if not (os.path.isdir(self.config.data_root_path)):
                raise Exception("Folder not found : %s" % self.config.data_root_path)
            if not (os.path.isdir(self.input_dir)):
                raise Exception("Folder not found : %s" % self.input_dir)
            copy_tree_with_absolute_symlinks(self.input_dir, dirname)

        # extract options
        configure_variant_options = self.variant["configure"]
        configure_case_option = f"--with-case={case_cpp_key}"
        configure_compiler_option = ""
        if tuning_flags != "":
            configure_compiler_option = f'FFLAGS="{tuning_flags} {tuning_flags_extra}"'

        # Add cppkeys
        configure_cppkeys_options = ""
        if "cppkeys" in self.case:
            reshape = []
            for key, value in self.case["cppkeys"].items():
                if value:
                    reshape.append(f"+{key}")
                else:
                    reshape.append(f"-{key}")
            if len(reshape) > 0:
                reshape_str = ",".join(reshape)
                configure_cppkeys_options = f"--with-keys={reshape_str}"

        # debug
        debug_option = ""
        if self.config.twin_chercker:
            debug_option += " --enable-twin-checker"

        # build command
        command = "%s/configure %s %s %s %s %s" % (
            croco_source_dir,
            configure_compiler_option,
            configure_case_option,
            configure_cppkeys_options,
            configure_variant_options,
            debug_option,
        )
        command = apply_vars_in_str(
            command, {"case": self.case, "tuning": self.config.host["tuning"]}
        )

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

        # copy log file in result dir
        if os.path.exists(os.path.join(self.dirname, "jobcomp.log")):
            shutil.copyfile(
                os.path.join(self.dirname, "jobcomp.log"),
                os.path.join(self.dirname_result, "jobcomp.log"),
            )

    def enable_rvtk_checking(self):
        # vars
        croco_build = self.croco_build
        variant_name = self.variant_name
        rvtk_ref_variant_name = self.config.variant_ref_name

        # enable what we need
        croco_build.cppdef_h_set_key("RVTK_DEBUG", True)

        # except for variant we enable RVTK_DEBUG_READ
        if variant_name == rvtk_ref_variant_name:
            croco_build.cppdef_h_set_key("RVTK_DEBUG_WRITE", True)
        else:
            croco_build.cppdef_h_set_key("RVTK_DEBUG_READ", True)

    def build(self, extra_info: str = "", force_rebuild: bool = False):
        if self.config.continue_on_error:
            try:
                self.build_internal(extra_info, force_rebuild)
            except Exception as e:
                Messaging.step_error("Fail to build !")
                self.config.report.report_status(
                    self.case_name, self.variant_name, "build", False, str(e)
                )
                return
        else:
            self.build_internal(extra_info, force_rebuild)
        self.config.report.report_status(
            self.case_name, self.variant_name, "build", True
        )

    def build_internal(self, extra_info: str = "", force_rebuild: bool = False):
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
        self.setup_variant()

        # some specific handling
        Messaging.step(f"Special command line configs...")
        if self.config.rvtk:
            self.enable_rvtk_checking()

        # compile
        self.compile()

    def symlink_all_rvtk_ref_files(self):
        # vars
        case_name = self.case_name
        rvtk_ref_variant_name = self.config.variant_ref_name

        # refdir
        ref_rundir = self.calc_rundir(rvtk_ref_variant_name, case_name)
        refdir = os.path.join(ref_rundir, "check_file_*")

        # debug_infos
        extra_files = [f"{ref_rundir}/debug_infos"]

        # link all
        for file in glob.glob(refdir) + extra_files:
            if not os.path.exists(os.path.basename(file)):
                os.symlink(file, os.path.basename(file))

    def run(self, extra_info: str = ""):
        # display
        Messaging.section(f"Running CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")

        # extract some vars
        command_prefix = self.variant["command_prefix"]
        environ = self.variant.get("environ", {}).copy()
        dirname = self.dirname
        runs = self.config.runs
        results = self.config.results
        host_tuning = self.config.host["tuning"]
        rvtk = self.config.rvtk
        is_rvtk_ref = self.variant_name == self.config.variant_ref_name

        # load tuning env & override if needed
        for key, value in host_tuning["environ"].items():
            environ[key] = value

        # apply vars
        command_prefix = apply_vars_in_str(
            command_prefix, {"case": self.case, "tuning": host_tuning}
        )

        # build end
        env_line = ""
        for var, value in environ.items():
            env_line += f'{var}="{value}" '

        # build command and run
        command = (
            "%s ../../../scripts/correct_end.sh %s ./croco TEST_CASES/croco.in.%s"
            % (env_line, command_prefix, self.case["case"].capitalize())
        )
        with move_in_dir(dirname):
            # link ref files
            if rvtk and not is_rvtk_ref:
                self.symlink_all_rvtk_ref_files()

            # run
            if self.config.continue_on_error:
                try:
                    result = run_hyperfine(
                        command,
                        logfilename="croco.log",
                        runs=runs,
                        verbose=self.config.verbose,
                    )
                except Exception as e:
                    Messaging.step_error("Error while running...")
                    result = {
                        "results": [
                            {
                                "command": command,
                                "mean": 0,
                                "median": 0,
                                "stddev": 0,
                                "min": 0,
                                "max": 0,
                                "times": [0],
                            }
                        ],
                        "error": str(e),
                    }
            else:
                result = run_hyperfine(
                    command,
                    logfilename="croco.log",
                    runs=runs,
                    verbose=self.config.verbose,
                )

        # status
        if "error" in result:
            self.config.report.report_status(
                self.case_name, self.variant_name, "run", False, result["error"]
            )
        else:
            self.config.report.report_status(
                self.case_name, self.variant_name, "run", True
            )

        # calc output
        os.makedirs(self.dirname_result, exist_ok=True)
        output_filename = f"{self.dirname_result}/result-{self.full_name}.json"

        # add some infos
        result["name"] = self.full_name
        result["original_path"] = output_filename
        result["hostname"] = platform.node()
        result["selected_host_config"] = self.config.use_host_config
        result["run_config"] = {
            "host": self.config.host,
            "case": self.case,
            "variant": self.variant,
        }

        # dump results
        with open(output_filename, "w+") as fp:
            json.dump(result, fp=fp, indent="\t")

        # copy log file in result dir
        for iteration in range(runs):
            logfilename = "croco_%i.log" % (iteration + 1)
            if os.path.exists(os.path.join(self.dirname, logfilename)):
                shutil.copyfile(
                    os.path.join(self.dirname, logfilename),
                    os.path.join(self.dirname_result, logfilename),
                )

    def apply_patches(self, patches):
        # vars
        croco_build = self.croco_build

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
                    patch_lines(file_filtered, [change])

    def setup_case(self):
        # apply the case paches
        case_name = self.case_name
        if self.variant["tuning_fflags"] == "debug":
            patches = self.case.get("patches_debug", {})
        else:
            patches = self.case.get("patches", {})
        cppkeys = self.case.get("cppkeys", {})

        # display
        Messaging.step(f"Apply case config : {case_name}")
        self.apply_patches(patches)

    def setup_variant(self):
        # apply the case paches
        variant_name = self.variant_name
        patches = self.variant.get("patches", {})
        cppkeys = self.variant.get("cppkeys", {})

        # display
        Messaging.step(f"Apply variant config : {variant_name}")
        self.apply_patches(patches)

    def check_one_file_from_seq_ref(self, filename: str) -> None:
        # extract vars
        case_name = self.case_name
        dirname = self.dirname
        variant_ref_name = self.config.variant_ref_name

        # if ref variant skip
        if self.variant_name == variant_ref_name:
            Messaging.step(
                f"Checking {case_name} / {filename} skiped for '{variant_ref_name}'"
            )
            return
        else:
            Messaging.step(f"Checking {case_name} / {filename}")

        # error
        seq_dir = self.calc_rundir(variant_ref_name, case_name)
        seq_file = os.path.join(seq_dir, filename)
        if not os.path.exists(seq_file):
            raise Exception(
                "Missing '%s', are you sure you ran case '%s' first to get a reference for checks ?"
                % (seq_file, variant_ref_name)
            )

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
            raise Exception(
                f"Missing '%s', are you sure you ran case %s in reference directory %s ?"
                % (ref_file, case_name, refdir)
            )

        # compare
        actual_file = f"{dirname}/{filename}"
        compare_netcdf_files(ref_file, actual_file)

    def check_one_file(self, filename: str):
        # compare
        if self.config.continue_on_error:
            try:
                self.check_one_file_internal(filename)
            except Exception as e:
                Messaging.step_error("Failed to compare !")
                self.config.report.report_status(
                    self.case_name, self.variant_name, "check", False, str(e)
                )
                return
        else:
            self.check_one_file_internal(filename)

        # report
        self.config.report.report_status(
            self.case_name, self.variant_name, "check", True
        )

    def check_one_file_internal(self, filename: str):
        # vars
        refdir = self.config.use_ref

        # swtich
        if refdir:
            self.check_one_file_from_ref_dir(filename)
        else:
            self.check_one_file_from_seq_ref(filename)

    def check(self):
        for filename in self.case["check_outputs"]:
            self.check_one_file(filename)

    def plotphy(self):
        dirname = self.dirname
        full_name = self.full_name
        filename = self.case["check_outputs"][0]  # to check
        actual_file = f"{dirname}/{filename}"
        result_dir = self.dirname_result

        if "plot_diag_script" in self.case:
            self.plot_diag_script = os.path.join(dirname, self.case["plot_diag_script"])
            Messaging.step(f"Plotting {full_name}")

            command = [
                self.plot_diag_script,
                "--no-show",
                "--makepng",
                "--file",
                actual_file,
                "--output-dir",
                result_dir,
            ]

            try:
                subprocess.run(command, check=True)  # Exécute la commande
                Messaging.step(
                    "Successfully executed %s with arguments --no-show --makepng"
                    % self.plot_diag_script
                )
                # report
                self.config.report.report_status(
                    self.case_name, self.variant_name, "plotphy", True
                )
            except subprocess.CalledProcessError as e:
                Messaging.step_error(
                    f"Error during execution of {self.plot_diag_script}: {e}"
                )
                self.config.report.report_status(
                    self.case_name, self.variant_name, "plotphy", False, str(e)
                )
        else:
            Messaging.step(f"No Plotting script provided")

    def dump_mesh(self, anim: bool = False):
        # vars
        dirname = self.dirname
        full_name = self.full_name
        result_dir = self.dirname_result

        # loop each
        for filename in self.case["check_outputs"]:
            actual_file = f"{dirname}/{filename}"
            dump_in_file = f"{result_dir}/mesh/"
            if anim:
                Messaging.step(f"Drawing mesh-anim {full_name} / {filename}")
            else:
                Messaging.step(f"Drawing mesh {full_name} / {filename}")
            mesh_drawer = MeshDrawer(actual_file, self.config, self.case_name)
            if anim:
                mesh_drawer.plot_all_variables_animation(dump_in_file)
            else:
                mesh_drawer.plot_all_variables(dump_in_file)

    def make_ref(self):
        # extract vars
        variant_ref_name = self.config.variant_ref_name
        refdir = self.config.build_ref
        case_name = self.case_name
        dirname = self.dirname
        patches_and_keys = {
            "case": {
                "patches": self.case.get("patches", {}),
                "keys": self.case.get("keys", {}),
            },
            "variant": {
                "patches": self.variant.get("patches", {}),
                "keys": self.variant.get("keys", {}),
            },
        }

        # check
        assert refdir

        # progress
        Messaging.step(f"Storing ref {case_name} into ref={refdir}")

        # we create only for seq
        if self.variant_name != variant_ref_name:
            return

        # create dir if not exist
        refdir_case = f"{refdir}/{case_name}"
        os.makedirs(refdir_case, exist_ok=True)

        # copy the files there
        for filename in self.case["check_outputs"]:
            orig_path = f"{dirname}/{filename}"
            dest_path = f"{refdir_case}/{filename}"
            Messaging.step(f"Storing file {case_name} / {filename}")
            shutil.copyfile(orig_path, dest_path)

        # also copy the config
        self.croco_build.copy_config(refdir_case, case_name, patches_and_keys)

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
            json.dump(self.case, fp, indent="\t")

        # dump git info
        with open(f"{refdir}/{case_name}/git.log", "w+") as fp:
            fp.write(
                subprocess.getoutput(
                    "git log HEAD~1..HEAD 2>/dev/null || echo 'No git !'"
                )
                + "\n"
            )

##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json

# termcolor
from termcolor import cprint


##########################################################
class Report:
    def __init__(self, config) -> None:
        self.config = config
        self.report = {}

    def report_status(
        self,
        case: str,
        variant: str,
        step: str,
        status: bool,
        message: str = "",
        log: str = "",
    ) -> None:
        full_name = f"{case}-{variant}"
        if not full_name in self.report:
            self.report[full_name] = {
                "case": case,
                "variant": variant,
                "status": {
                    "build": {"status": None},
                    "run": {"status": None},
                    "check": {"status": None},
                    "plotphy": {"status": None},
                },
            }
        self.report[full_name]["status"][step] = {
            "status": status,
            "message": message,
            "log": log,
        }

    def save(self) -> None:
        status_report_file = os.path.join(self.config.results, "status_report.json")
        with open(status_report_file, "w+") as fp:
            json.dump(self.report, fp, indent="\t")

    # To detect if report contains False values
    def contains_false(self):
        for key_case, value_case in self.report.items():
            status_dict = value_case.get("status")
            for key, value in status_dict.items():
                if isinstance(value, dict) and value.get("status") is False:
                    return True
            return False

    def get_color(self, status):
        if status is True:
            return "green"
        elif status is False:
            return "red"
        elif status is None:
            return "cyan"
        else:
            return "white"

    def print(self) -> None:
        print(
            "------------------------------- REPORT SUMMARY ------------------------------------"
        )

        # header
        full_name = "CASE-VARIANT"
        status_build = "Build"
        status_run = "Run"
        status_check = "Check"
        status_plotphy = "Plotphy"
        print(
            "%s %s %s %s %s"
            % (
                f"{full_name:<40}",
                f"{str(status_build):<8}",
                f"{str(status_run):<8}",
                f"{str(status_check)}",
                f"{str(status_plotphy)}",
            )
        )
        print(
            "..................................................................................."
        )

        for full_name, infos in self.report.items():
            status_build = infos["status"]["build"]["status"]
            status_run = infos["status"]["run"]["status"]
            status_check = infos["status"]["check"]["status"]
            status_plotphy = infos["status"]["plotphy"]["status"]

            cprint(f"{full_name:<40} ", "white", end="")
            cprint(f"{str(status_build):<8} ", self.get_color(status_build), end="")
            cprint(f"{str(status_run):<8} ", self.get_color(status_run), end="")
            cprint(f"{str(status_check):<8} ", self.get_color(status_check), end="")
            cprint(f"{str(status_plotphy):<8}", self.get_color(status_plotphy))

        print(
            "-----------------------------------------------------------------------------------"
        )

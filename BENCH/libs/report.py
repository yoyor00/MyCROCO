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

    def report_status(self, case: str, variant: str, step: str, status: bool, message: str = '', log: str = '') -> None:
        full_name = f"{case}-{variant}"
        if not full_name in self.report:
            self.report[full_name] = {'case': case, 'variant': variant, 'status': {
                'build': {'status': None},
                'run': {'status': None},
                'check': {'status': None}
            }}
        self.report[full_name]['status'][step] = {
            'status': status,
            'message': message, 
            'log': log
        }

    def save(self) -> None:
        status_report_file = os.path.join(self.config.results, "status_report.json")
        with open(status_report_file, "w+") as fp:
            json.dump(self.report, fp, indent='\t')

    def print(self) -> None:
        print("------------------------------- REPORT SUMMARY ------------------------------------")

        # header
        full_name = 'CASE-VARIANT'
        status_build = 'Build'
        status_run = 'Run'
        status_check = 'Check'
        print(f"{full_name:<40}    {str(status_build):<8}  {str(status_run):<8}      {str(status_check)}")
        print("...................................................................................")

        for full_name, infos in self.report.items():
            status_build = infos['status']['build']['status']
            status_run = infos['status']['run']['status']
            status_check = infos['status']['check']['status']
            if status_build and status_run and status_check:
                color = 'green'
            elif not status_build or not status_run or not status_check:
                color = 'red'
            else:
                color = 'orange'
            cprint(f"{full_name:<40}    {str(status_build):<8}  {str(status_run):<8}      {str(status_check)}", color)

        print("-----------------------------------------------------------------------------------")
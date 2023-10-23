##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import platform
import subprocess
import json
from datetime import datetime

##########################################################
def run_and_get_first_line(command: str) -> str:
    try:
        return subprocess.check_output("mpirun --version").split('\n')[0]
    except:
        return "UNAVAILABLE"

##########################################################
def get_lscpu_as_json() -> dict:
    data = subprocess.check_output(['lscpu', '-J'])
    return json.loads(data)

##########################################################
def gen_system_info() -> dict:
    # dict
    infos = {}

    # date
    infos['date'] = datetime.now().strftime("%d/%m/%y")
    infos['time'] = datetime.now().strftime("%H:%M:%S")

    # get hostname
    system = infos['plateforme'] = {}
    system['hostname'] = platform.node()
    system['arch'] = platform.architecture()
    system['system'] = platform.system()
    system['uname'] = platform.uname()
    
    # lscpu
    lscpu = get_lscpu_as_json()['lscpu']
    infos['lscpu'] = lscpu
    
    for entry in lscpu:
        if 'children' in entry:
            system['processor'] = entry['children'][0]['data']

    # get mpi type
    soft = infos['soft'] = {}
    soft['mpi'] = run_and_get_first_line("mpirun --version")
    soft['gcc'] = run_and_get_first_line("gcc --version")
    soft['gfortran'] = run_and_get_first_line("gfortran --version")
    soft['nvfotran'] = run_and_get_first_line("nvfortran --version")
    soft['cmake'] = run_and_get_first_line("cmake --version")

    # ok
    return infos

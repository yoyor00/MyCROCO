##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
class Messaging:
    '''Define the messaging functions to print progress.'''

    @staticmethod
    def section(message: str) -> None:
        '''Print progress with section'''
        print(f" - {message}")

    @staticmethod
    def command(message: str) -> None:
        '''Before executing a command'''
        print(f"     >>>> {message}")

    @staticmethod
    def step(message: str) -> None:
        '''Print progress with a step'''
        print(f"     ---- {message}")

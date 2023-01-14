######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
cmake_minimum_required(VERSION 3.1)

######################################################
project(CROCO
        VERSION 1.3.0
        DESCRIPTION "Oceanic simulation"
        LANGUAGES Fortran C)

######################################################
add_subdirectory(OCEAN)

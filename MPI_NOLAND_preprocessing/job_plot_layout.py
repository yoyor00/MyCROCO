#!/bin/bash
set -x
cd ${BRIDGE_MSUB_PWD}
pwd
#==============================================================================
#                      header to submit a MPI batch job 
#===============================================================================

# Type ccc_mprun -h for an updated and complete documentation.

#MSUB -r giga3       # request name
#MSUB -n 1         # Total number of mpi task to use # for GIGATL3_48x48 
#MSUB -T 600         # elapsed time limit in seconds

#MSUB -o giga3%I.txt
#MSUB -e giga3%I.txt

#MSUB -A gen7638
#MSUB -q xlarge 
#MSUB -@ sebastien.theetten@ifremer.fr:begin,end       # only for thin nodes 
#===============================================================================

module load python

#ccc_mprun python plot_layout.py gigatl6_grd.nc GIGATL6-037x027_0708
ccc_mprun python plot_layout.py gigatl1_grd.nc GIGATL1-0100x0100_6589


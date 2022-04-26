Tutorial
========

.. _tutorial:

Setup on datarmor
-----------------

You can download the latest version of CROCO with BIOLink in the folder :

/home1/homework/gkoenig/croco/

This folder can be downloaded in your working directory with the command :

cp -r /home1/homework/gkoenig/croco/ .

Exploration of the folder 
-------------------------

In this folder you will find the usual folders of CROCO, with three supplementary folders :
BIOLink, BLOOM and VILAINE. Those folders contain the coupler, the BLOOM NPZD model and the the VILAINE configuration for CROCO and BLOOM.

The codes in themselves are not to be modified, but feel free to explore them. Configuration files are stored in the VILAINE folder to modify the parameters of the BLOOM model. There are seven of those files.

.. list-table:: Parameter files to be modified for using BIOLink
   :widths: 25 25 50
   :header-rows: 1

   * - File
     - Position
     - Model
   * - cppdefs.h
     - VILAINE/
     - CROCO
   * - param.h
     - VILAINE/
     - CROCO
   * - croco.in
     - VILAINE/
     - CROCO
   * - parasubstance_BLOOM.txt
     - VILAINE/MUSTANG_NAMELIST
     - SUBSTANCE
   * - parabloom.txt
     - VILAINE/FIC_NAMELIST
     - BLOOM
   * - param_BLOOM.h
     - BLOOM/
     - BLOOM
   * - vardiag.dat_debase
     - BLOOM/
     - BLOOM



The different files are dedicated to activate different options in the models. The files of CROCO change the number of tracers, the size of the grid and the precompilation (cpp) keys.

The precompilation keys are grouped in cppdefs after the define SUBSTANCE and define BIOLink part :

.. figure:: /cppdefsBIOLink.png
   
   Example of the modifications of the cppkeys for BIOLink
   
The SUBSTANCE and TRACERS keys are required because handles the biological/sediment variables as tracers. The SUBSTANCE submodel is required for the variables that can change the density of seawater or for particulate variables that can sink.
  
The file param.h also have to be partly modified.
  
.. figure:: /paramh.png
   
   Example of modification of param.h for BIOLink
   
Here we declare three new variables to count the number of tracers: ntrc_sedim, ntrc_biol and ntrc_fix which are the number of variables of the sediments, biological and fixed. The ntrc_biol and ntrc_fix numbers are determined in the file param_bloom.h

The croco.in file contains the post-compilation parameters. This includes the position of forcing files, variables recorded in the history and averaged files and the diffusion coefficients values. The number of tracers as to be modified in several places to correspond to the number declared with BIOLink.

The first changes appear around the line 43.

primary_history_fields: zeta UBAR VBAR U V wrtT(1:NT)
                         T    T    T   T T  30*T

and 50 :

primary_averages: zeta UBAR VBAR U V wrtT(1:NT)
                   F    F    F   F F  30*F

The lines after "primary_history_fields" and "primary_averages" are not read. The modifications are to be made under the wrtT(1:NT) where the numbers 30*T and 30*F are to be changed to the number of tracers ( including the salinity and temperature). If it is not known exactly, any number superior to the number of tracers will work.

The same modifications are to be made for the vertical and horizontal diffusion coefficients.

tracer_diff2: TNU2(1:NT)
               30*0.d-2 [m^2/sec for all]
tracer_diff4: TNU4(1:NT)
               30*0.d11 [m^4/sec for all]
vertical_mixing: Akv_bak, Akt_bak [m^2/sec]
                   0.d0    30*0.d0


In the same way, the first number of values written under TNU2, TNU4 and Akt_bak are to be adjusted to the number of tracers. 

Next one needs to change the path to the parameter file of BLOOM, which will be described further.

substance: input file (parasubstance)
MUSTANG_NAMELIST/parasubstance_BLOOM.txt

Finally, the number of tracers as to be modified one last time to corresponds to the river inputs. Those are in the lines 115 to 119 :

psource_ncfile: Nsrc Isrc Jsrc Dsrc qbardir Lsrc Tsrc
runoff file name
CROCO_FILES/croco_runoff_vilaine_bloom_var.nc
                 2    
                       167 56 0 -1   30*T 20.0 15.0 Loire
                       91 99 0 -1    30*T 20.0 15.0 Vilaine_arzal
                       

Once again, the number of tracers before the *T have to be modified to agree with the number of tracers.

Next we have to modify the file parasubstance_BLOOM.txt for substance. THe substance submodel is a code developed to interface the sedimentary model MUSTANG with the hydrodynamical model CROCO. It has been extended to handle non-sedimentary variables. It is used in BIOLink to defined characteristics of various tracer variables. Those characteristics are defined in the file parasubstance_BLOOM.txt. It contains a set of namelists which describe the tracers. The first namelist contains the number of variables and must be coherent with the number declared in croco.in and param_bloom.h


&nmlnbvar
nv_dis=4 ! number of dissolved susbtances
nv_dis=4 ! number of Non Constitutive Particulate subtances
nv_ncp=9 ! number of fixed susbtances (not advected)
nv_fix=3 ! number of of particulate susbtances sorbed on an other particule
nv_sorb=1 ! number of of benthic susbtances
nv_bent=0 /

The different kind of variables correspond to different behaviors inherited from the sedimentary model. For example if the variables sink or not.

For BLOOM there are three configuration files. The first, parabloom.txt contains parameters for the model BLOOM itself. The second, param_BLOOM.txt contains the number of variables used as tracer in BLOOM and finally, vardiag.dat_debase contains the diagnostic variables of BLOOM. Those files have not been modified for the model in VILAINE. 

Compilation
-----------

The Makefile has already been modified and everything for the compilation on datarmor has been included in the script batch_comp_datarmor. You can launch it with the command :

qsub batch_comp_datarmor

First run
---------

The model can then be run with the datarmor command 

qsub batch_mpt

Analysis of outputs
-------------------

The outputs will appear in a history file given in the croco.in file.

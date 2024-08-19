DJL matlab package

To run Croco SOLITON_DJL, you will have to run make_djl.m Matlab script.
This script is using a couple of routines from the croco_tools library. Do set
croco_tools_path to ensure the correct behaviour of make_djl.m.

crocotools_param_djl.m provides the parameters for the initialization of SOLITON_DJL 
Croco configuration.


You have to set a few parameters from the Matlab file to adapt to your domain size/depth. 
H  : depth of your domain
L  : length of your domain
NX : number of x grid points (+2)
NZ : number of vertical levels (Z-grid)
f_shift : the distance along x-direction for shifting the soliton pattern 
          (caution when domain is refined!)

Note that an increase of resolution is required later in the code, there you should provide 
the real NX, NZ for your croco grid 

The DJL algorithm asks for an APE for the wave as an initilization. 

The DJL algorithm asks also for a seawater equation of state, identical to yours.

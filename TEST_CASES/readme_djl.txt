DJL matlab package

To provide Croco SOLITON_DJL run input parameters, you wil have to run make_djl.m Matlab
code.
This code is using a couple of routines from the croco_tools library. Do set
croco_tools_path to ensure the corect behaviour of make_djl.m.

The crocotools_param_djl.m provide the parameters for the initialisation of the SOLITON_DJL 
croco configuration.


You have to set few parameter from the matlab file to adap to your domain size/depth. 
H: the depth of you domain
L:The lenght of your domain
NX : number of x grid point
NZ : number of vertical level (Z-grid)
f_shift : the distance along x-direction for shifting the soliton pattern 

Note that an increase of resolution is asked later in the code, there you shoul provide 
the good NX, NZ for you croco grid 

The DJL algorithm ask an APE for the wave as initilisation. 

The DJL algorithm ask also a seawater equation of state, identical to your croco.in one  

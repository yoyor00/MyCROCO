CROCO v1.1
-----------
Released date : 25 October 2019

Previous release : CROCO v1.0 (June 2018)


Reminders:
==========
CROCO sources and CROCO_TOOLS (the follow-on of ROMS_TOOLS) are now distributed separately (for croco_tools releases, see associated tab at  https://www.croco-ocean.org/download/croco-project/ ). ROMS_AGRIF  is not maintained anymore and we strongly encourage ROMS_AGRIF users to switch  to CROCO. CROCO version available directly from the git repository is a unstable development  version. Standard users should use the stable one downloaded from the web site. 
 
New in v1.1 :
=============

- architecture : slight changes
	- CROCO doesn't come with a default Run directory anymore, first step as a user will be to edit the script create_run.bash to create your own
	- as a consequence default TEST_CASES are now stored in croco/TEST_CASES and the scipts for running under croco/SCRIPTS

- parallelisation : in MPI (not OPENMP), CROCO has now the capacity to avoid compuation on land only processors. A pre-processing step is required. The corresponding tool is located under croco/MPP_PREP : the namelist has to be filled ith the maximum number of CPUs available, then compile and execute the code. It returns the optimal parameters for NP_XI, NP_ETA and NNODES.
# define MPI
# define MPI_NOLAND   (default undef)

- outputs : NetCDF4 parallel cpabilities are available. When running with MPI activated, writing in a single file (without PARALLEL_FILES) is faster
# define MPI
# define NC4PAR   (default undef)
            XIOS has been updated to XIOS2.5. A sample of the new xml files are available in croco/XIOS. There is no backward compatibility with XIOS1. See https://forge.ipsl.jussieu.fr/ioserver/wiki/documentation 
# define MPI
# define XIOS   (default undef)

- diagnostics : several trends computation added. Each of them has its corresponding section in croco.in.
                vorticity : 
# define  DIAGNOSTICS_VRT
                eddy kinetic energy :
# define  DIAGNOSTICS_EK
                potential vorticity  
# define  DIAGNOSTICS_PV

- grid : possibility of refining locally the resolution of meshes to take into account section reduction (example : schematic river canal in a realistic configuration : keep realistic sections based on a 500m meshgrid size)
# define REDUC_SECTION  (default undef)

- bottom friction : Spatialisation of Z0 coefficient
# define Z0B_VAR   (default undef)
					Use sub time steps for computation of botoom friction. In this case, bottom friction is applied in step3d_fast routine
# default BSTRESS_FAST  (default depends on the configurations) 

- functional : Use of calendar to define start and end date : The user can choose the beginning and the end of the simulation in croco.in , then the number of time steps will be calculated. Also print the date at each time step in screen log.
# define USE_CALENDAR   (default undef)
with in croco.in the additional lines
  run_start_date:
    01/04/2015 00:00:00
  run_end_date:
    12/04/2015 00:00:00

- atmospheric forcing : Manage MeteoFrance inputs (Arome/Arpege) for online atmospheric reading and interpolation
# define AROME        (default undef)
						Read atmospherical pressure from meteo file and use it in bulk flux and take into account atmospherical pressure gradient in the equations
# define READ_PATM	  (default undef)
						Current feedback added in bulk and forcing cases. Three  possibilities: feedback on the stress (CFB_STRESS, default) using a coefficient based on the wind speed, or based on the wind stress (CFB_STRESS2), or feedback on the wind speed (CFB_WIND, less physical)
# define SMFLUX_CFB (default undef)

- ocean/wave coupling : add wave average outputs. Add wave-to-ocean stress accounting for wave model bulk

- tides : Introduction of a new method for harmonic composition : use of Simon (SHOM) method to build tide elevation from harmonic constituents.
# define TIDES_MAS    (default undef)

- open boudaries : reduced-form equation for barotropic velocities at boundary conditions in case of velocity tides are not available (not defined UV_TIDES). It's also useful for test cases to only force the model with an analytical sea surface elevation (for example a M2 tide). In both cases it replaces M2_FRC_BRY (and this cpp key should be set to False to avoid confusion) and it works with M2FLATHER and M2CHARACT obc conditions. This development comes from Roms-Rutgers source code.
This option has been tested in a tidal flat test case for each open boundary condition (north, south, west and east)
# define OBC_REDUCED_PHYSICS (default undef)
					OBC_M2FLATHER has been suppressed. Use OBS_M2CHARACT instead

- vertical mixing : GLS scheme has been rewritten. GLS_MIX2017 cpp key as been suppress
#define GLS  (default depends on the configurations, if activated KEPSILON model and  CANUTO stability function by default)
Dedicated configuration : SINGLE_COLUMN, a 1D vertical model with several subsettings to test.

- test cases : addtion or renaming
I_SOLITON : was previously GRAV_ADJ with GRAV_ADJ_SOLITON activated
FLUME : now SANDBAR
TS_HADV_TEST : new case to test horizontal advection for tracers



 _   _                      __                 _
| | | | __ ___   _____     / _|_   _ _ __     | |
| |_| |/ _` \ \ / / _ \   | |_| | | | '_ \    | |
|  _  | (_| |\ V /  __/   |  _| |_| | | | |   |_|
|_| |_|\__,_| \_/ \___|   |_|  \__,_|_| |_|   (_)

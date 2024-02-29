#-------------------------------------------------------------------------------
# OCE
#-------------------------------------------------------------------------------
# namelist [Info: grid size is directly read in oce_compile.sh and cpl_nam.sh ]

# Online Compilation
#-------------------
#!!!!!!! IMPORTANT NOTE !!!!!!!
# If activated, creates croco executable depending on this namelist. 
#   - In param.h it modifies the grid size, the number of procs in x and y direction with those given in myjob.sh
#   - In cppdefs.h it modifies the following options with informations given below
#       MPI, OA_COUPLING, OW_COUPLING, MRL_WCI, 
#       XIOS, LOGFILE, MPI_NOLAND,
#       AGRIF, AGRIF_2WAY, 
#       BULK_FLUX, ONLINE, AROME, ARPEGE, ERA_ECMWF
#       FRC_BRY, CLIMATOLOGY
#       TIDES, PSOURCE, PSOURCE_NCFILE, PSOURCE_NCFILE_TS
# Other changes of parameterizations, numerical schemes, etc should be made "by hand" in CROCO_IN/cppdefs.h.base
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export ONLINE_COMP=0

# Time steps
#-----------
export DT_OCE=3600
export NDTFAST=60

# domains
#--------
export AGRIFZ=0
export AGRIF_2WAY="FALSE"

# MPI_NOLAND
#-----------
export MPI_NOLAND="FALSE"
export MY_NODES=18 # ONLY if MPI_NOLAND is "TRUE". It replaces NP_OCEX*NP_OCEY

# wave coupling
#--------------
export OW_COUPLING_FULL="FALSE" # use all terms from the wave model (stokes, tke flux, bernouilli head, orbital velocities) if TRUE or use their monochromatic approximation if FALSE
export WAVE_SMFLUX="FALSE" # use stress from the wave model if TRUE or that from the bulk formulation if FALSE. Not used when cpl with atm

# forcings
#---------
export ini_ext='ini_SODA' # ini extension file (ini_SODA,...)
export bdy_ext='bry_SODA' # bry extension file (clm_SODA,bry_SODA,...)

export surfrc_flag="TRUE" # Flag if surface forcing is needed (FALSE if coupling with the atmosphere, TRUE otherwise)
export interponline=0 # switch (1=on, 0=off) for online surface interpolation. Only works with MONTHLY input files!
export frc_ext='blk_ERA5' # surface forcing extension(blk_ERA5, frc_ERA5,...). If interponline=1 precise the type (ERA_ECMWF or AROME,  [CFSR by default], names as cppkey name in croco)

export tide_flag="FALSE" # the forcing extension must be blk_??? otherwise tide forcing overwrites it 
export river_flag="FALSE"

# output settings
#----------------
#!!!!!!! WARNING !!!!!!!
# When XIOS is activated the following values (for the model) are not taken into account
#!!!!!!!!!!!!!!!!!!!!!!!
export oce_his_sec=86400     # history output interval (in number of second) 
export oce_avg_sec=86400     # average output interval (in number of second) 



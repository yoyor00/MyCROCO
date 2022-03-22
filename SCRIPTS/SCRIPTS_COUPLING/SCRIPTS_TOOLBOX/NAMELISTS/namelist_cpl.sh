#-------------------------------------------------------------------------------
# CPL
#-------------------------------------------------------------------------------
# namelist
#!!!! WARNING !!! 
# If you are using pre-built weight. Please modified manually the netcdf used in the namcouple
# add .smtho2a if you want to smooth
export namcouplename=namcouple.base.${RUNtype} 

# coupling frequency
export CPL_FREQ=3600

# Start from existing condition
export CPL_restart="FALSE" # If TRUE: initialize CPL field with an history field
 
# Files to create oasis restart (leave empty if none). File need to be in *FILES_DIR. Default 1st time step
export oce_rst_file=""   
export atm_rst_file=""
export wav_rst_file""

# Use pre-build weight files (need to be in OASIS_FILES)
#!!!! WARNING !!!!
# At the moment only workd with 1-atm and 1-oce domain
#!!!!!!!!!!!!!!!!!
export WEIGHT_FLAG=0 
export weight_o2a=""
export weight_a2o=""


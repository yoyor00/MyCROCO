#-------------------------------------------------------------------------------
# CPL
#-------------------------------------------------------------------------------
# namelist
export namcouplename=namcouple.base.${RUNtype}

# coupling frequency
export CPL_FREQ=3600

# Start from existing condition
export CPL_restart="FALSE" # If TRUE: initialize CPL field with an history field
 
# Files to create oasis restart (leave empty if none). File need to be in *FILES_DIR. Default 1st time step
export oce_rst_file=""   
export atm_rst_file=""
export wav_rst_file""


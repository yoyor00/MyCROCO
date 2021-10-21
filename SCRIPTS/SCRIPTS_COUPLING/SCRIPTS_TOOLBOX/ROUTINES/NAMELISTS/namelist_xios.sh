#-------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------
# namelist
export FILIN_XIOS="iodef.xml context_croco.xml file_def_croco.xml domain_def_croco.xml field_def_croco.xml"
#export FILIN_XIOS="iodef.xml context_croco.xml file_def_croco.xml domain_def_croco.xml field_def_croco.xml context_wrf.xml file_def_wrf.xml" # files needed for xios. Need to be in INPUTDIRX (cf header.*)

# files
export ATM_XIOS_NAME="wrfout_inst" # All the names you have defined in your xml file
export OCE_OUTPUT_PREFIX="croco"
export OCE_XIOS_NAME="${OCE_OUTPUT_PREFIX}_1d_aver ${OCE_OUTPUT_PREFIX}_1h_inst_surf ${OCE_OUTPUT_PREFIX}_6h_his"


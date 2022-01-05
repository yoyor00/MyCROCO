#-------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------

export ONLINE_XML="FALSE" # Process xml online

# name of xml files (defined in file_def_*). Here are default names in xml files.
export ATM_XIOS_NAME="wrfout wrf3d_1D wrf3d_1H" # All the names you have defined in your xml file
export OCE_OUTPUT_PREFIX="croco"
export OCE_XIOS_NAME="${OCE_OUTPUT_PREFIX}_1d_aver ${OCE_OUTPUT_PREFIX}_1h_inst_surf ${OCE_OUTPUT_PREFIX}_6h_his"


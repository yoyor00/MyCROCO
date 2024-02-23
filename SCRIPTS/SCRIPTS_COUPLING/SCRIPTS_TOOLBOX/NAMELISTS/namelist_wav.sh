#-------------------------------------------------------------------------------
# WAV
#-------------------------------------------------------------------------------

# Time steps
#-----------
export DT_WAV=3600     # TMAX = 3*TCFL
export DT_WW_PRO=1200  # TCFL = 0.8 x dx/(g/fmin4pi) with fmin=0.0373 => 3-4 % of dx
export DT_WW_REF=1800  # TMAX / 2
export DT_WW_SRC=10    # TSRC = usually 10s  (could be between 5s and 60s)

# Grid size
#----------
export wavnx=41 ; export wavny=42

# Parameter
#---------- 
export hmin=75; # e.g. minimum water depth in CROCO (will be used to delimit coastline in WW3)

# Forcing files
#--------------
# forcing file(s) PREFIX list (leave empty if none)
# input filenames are supposed to be in the form: PREFIX_Y????M??.nc
export forcin=() 
# name of ww3_prnc.inp extension, e.g wind, see in WW3_IN directory
export forcww3=() 

# Boundary files
#---------------
# prefix for boundary files (leave empty is none)
export bouncin= 

# Output settings
#----------------
export wav_int=21600  # output interval (s)
export flagout="TRUE" # Keep (TRUE) or not (FALSE) ww3 full output binary file (out_grd.ww3)


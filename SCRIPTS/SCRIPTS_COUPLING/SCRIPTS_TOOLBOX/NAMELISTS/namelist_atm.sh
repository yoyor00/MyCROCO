#-------------------------------------------------------------------------------
# ATM
#-------------------------------------------------------------------------------
# namelist
export atmnamelist=namelist.input.base.complete

# Time steps
export DT_ATM=100 #100   # 100 90 75 72 60 45

# Grid size
#[ Grid size should be already put in the namelist. When coupled it is directly read in cpl_nam.sh ]

# domains
export NB_dom=1 # Number of coupled domains
export wrfcpldom='d01'
export nestfeedback="TRUE" # 1 way (FALSE) or 2 Way (TRUE) nesting
export onlinecplmask="TRUE" # Build default mask (depending on the nb of atm and oce domains)
# Boundaries interval 
export interval_seconds=21600 # interval ( in sec ) of the latteral input
export auxinput4_interval=360 # interval ( in min ) of bottom input
export ptop=5000
export nbmetsoil=4
export nbmetlevel=38
# fdda options
switch_fdda=0
nudge_coef=0.0003
nudge_interval_m=360
nudge_end_h=144

# output settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                          WARNING                                       ! 
# When XIOS is activated the following values (for the model) are not taken into account !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export atm_his_h=6                        # output interval (h)
export atm_his_frames=1000 # $((31*24))          # nb of outputs per file
export atm_diag_int_m=$((${atm_his_h}*60))  # diag output interval (m)
export atm_diag_frames=1000     # nb of diag outputs per file


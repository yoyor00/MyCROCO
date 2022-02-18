#-------------------------------------------------------------------------------
# ATM
#-------------------------------------------------------------------------------
# What kind of run
export ATM_CASE="DEFAULT"  # DEFAULT or MOVING_NEST

# namelist
export atmnamelist=namelist.input.base.complete

# Time steps
export DT_ATM=100

# Grid size
#[ Grid size should be already put in the namelist. When coupled it is directly read in cpl_nam.sh ]

# domains
export NB_dom=1 # Number of coupled domains
export wrfcpldom='d01'
export nestfeedback="TRUE" # 1 way (FALSE) or 2 Way (TRUE) nesting
export onlinecplmask="TRUE" # Erase existing CPLMASK and build default mask (depending on the nb of atm and oce domains)
# Boundaries interval 
export interval_seconds=21600 # interval ( in sec ) of the latteral input
export auxinput4_interval=360 # interval ( in min ) of bottom input
export nbmetsoil=4
export nbmetlevel=38
# fdda options
export switch_fdda=0 # To activate fdda
export nudgedom="1" # which kind of nudging (0,1,2) see WRF_IN/README.namelist. If parent + nest "1 1"
export nudge_coef="0.0003" # nudge coef. Need to be the same size than nudge
export nudge_interval_m="360" # time interval (in min) between analysis times 
export nudge_end_h="144" # time (in hours) to stop nudging after start of forecast
# physics
export isftcflx=0 # Cd formulation for tropical storm application (default 0, wave cpl =5)

##### MOVING NEST CASE ####
export num_mv_nest=1 # number of moving nest

# Nest informations (if several nest, the following variables need to have the format "1st_nest 2nd_nest" ) 
export ref_coef="3" # refinement coef for nest
export ew_size="283" # nest size in east-west dim ([multiple of ref_coef] + 1)
export ns_size="295" # nest size in north-south dim ([multiple of ref_coef] + 1)
export i_prt_strt="580" # where nest is starting in parent's grid x-dim 
export j_prt_strt="59" # where nest is starting in parent's grid y-dim

# Tracking parameter
vortex_interval=5 # When to update vortex position
max_vortex_speed=40 # Used to compute the search radius for the new vortex center position
corral_dist=8 # The closest distance between child and parend boundary (in parent grid cell)
track_level=50000 # The pressure level (in Pa) where the vortex is tracked
time_to_move=0 # The time (in minutes) until the nest is moved (at the beginning)

# output settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                          WARNING                                       ! 
# When XIOS is activated the following values (for the model) are not taken into account !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export atm_his_h=6                        # output interval (h)
export atm_his_frames=1000 # $((31*24))          # nb of outputs per file
export atm_diag_int_m=$((${atm_his_h}*60))  # diag output interval (m)
export atm_diag_frames=1000     # nb of diag outputs per file


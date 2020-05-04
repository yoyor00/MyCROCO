#!/bin/bash

#set -x
#-----------------------
source CONFIGURE_GLOBAL
#-----------------------

#-1 Set-up the CVTK_FAST environement for reg, ana and vort type
./Scripts_reg/create_link_master_reg.sh
./Scripts_ana/create_link_master_ana.sh
./Scripts_vort/create_link_master_vort.sh

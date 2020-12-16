G. Cambon December 2020, IRD

Warning
=======
- TAKE CARE : cppdefs.h_exactrst_debugwrite and cppdefs.h_exactrst_debugread only differ by the
cppkeys

for cppdefs.h_exactrst_debugwrite => # define RVTK_DEBUG_READ
for cppdefs.h_exactrst_debugread  => # define XXXXXRVTK_DEBUG_READ

- TAKE CARE:  the cppdefs.h_* are coherent with croco/OCEAN/cppdefs.h

- TAKE CARE: in croco.in.rst_seg12_PR, you must have
restart:          NRST, NRPFRST / filename
                   6      2
    OUT_seg12/croco_rst.nc

- TAKE CARE in croco.in,  you must have
initial: NRREC  filename
          2
    CROCO_FILES/croco_rst.00000.nc


Procedure
=========
Using the configuration BENGUELA_VHR [ ftp://ftp.ifremer.fr/ifremer/croco/ForCI-INRIA/VHR_CROCO_FILES_BCK ]

./compile_croco_PR.bash
./test_debug_perfrst.bash

If using another configuration, adapt the input file paths in the croco.in*


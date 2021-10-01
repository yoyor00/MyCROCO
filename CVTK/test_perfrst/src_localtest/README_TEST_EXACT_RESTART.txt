G. Cambon December 2020, IRD
----------------------------

The test procedure of the EXACT_RESTART make use of the cpp-key RVTK_DEBUG_PERFRST

1- we do run restarted at it=6 ( run1; it=1-6 then restart it=6-12). We wrote in a binary file some critical arrays between it=6 and 12.  keys needed : # define RVTK_DEBUG_READ + # define RVTK_DEBUG_PERFRST + # define EXACT_RESTART

2- we do a continuous run (it=1-12) run and we compare read and compare  the array from step 6 to 12
keys needed  # undef RVTK_DEBUG_READ  (/ or # define XXXXXRVTK_DEBUG_READ) + # define + #define RVTK_DEBUG_PERFRST + # define EXACT_RESTART

For more info, see debug.F

Procedure
=========
Using the configuration BENGUELA_VHR [ ftp://ftp.ifremer.fr/ifremer/croco/ForCI-INRIA/VHR_CROCO_FILES_BCK ]

./compile_croco_PR.bash
./test_debug_perfrst.bash

If using another configuration, adapt the input file paths in the croco.in*

Warning
=======
- TAKE CARE : cppdefs.h_exactrst_debugwrite and cppdefs.h_exactrst_debugread only differ by the
cppkeys

- for cppdefs.h_exactrst_debugwrite => # define RVTK_DEBUG_READ
- for cppdefs.h_exactrst_debugread  => # define XXXXXRVTK_DEBUG_READ
- for both, it needs 
      => # define RVTK_DEBUG_PERFRST
      => # define EXACT_RESTART

- TAKE CARE:  the cppdefs.h_* are coherent with croco/OCEAN/cppdefs.h

- TAKE CARE: in croco.in.rst_seg12_PR, you must have

restart:          NRST, NRPFRST / filename
                   6      2
    OUT_seg12/croco_rst.nc

- TAKE CARE in croco.in,  you must have

initial: NRREC  filename
          2
    CROCO_FILES/croco_rst.00000.nc



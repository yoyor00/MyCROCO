#!/bin/bash

#set -x

rm -Rf croco_EXACT_RST*
rm -Rf Compile_read Compile_write

ln -sf cppdefs.h_exactrst_debugwrite cppdefs.h
#ln -sf param.h_debugwrite param.h
./jobcomp
mv croco croco_EXACT_RST_DBWRITE

cp -Rf Compile Compile_write

ln -sf cppdefs.h_exactrst_debugread cppdefs.h
#ln -sf param.h_debugread param.h
./jobcomp
mv croco croco_EXACT_RST_DBREAD
cp -Rf Compile Compile_read

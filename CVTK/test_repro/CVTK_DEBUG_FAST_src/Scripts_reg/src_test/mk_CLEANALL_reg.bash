#!/bin/bash

for i in `ls Configure_Test_reg/*` ; do
    ii=`echo $i | cut -d/ -f2-`
    echo '=> Clean the dir: '$ii
    rm -Rf $ii
done

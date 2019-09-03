#!/bin/bash

#set -e
set -u
#set -x

source CONFIGURE_GLOBAL

./git_process.bash

for testconf in `ls -1 ./Configure_Test/ `;do
  [ -d $testconf ] && rm -rf 	$testconf 
  echo -e ${FMT_BLUEBLD}"=============================="${FMT_ORD}
  echo -e ${FMT_BLUEBLD}"TESTING $testconf :"${FMT_ORD}
  echo -e ${FMT_BLUEBLD}"=============================="${FMT_ORD}
  rm -rf $testconf
  ./mk_TestDIR_ana.bash $testconf
  echo "  "  
done

found=0
for testconf in `ls -1 ./Configure_Test/ `; do
  mytest=$(ls $testconf/jobcomp_OPENMP* > /dev/null 2>&1)
  if [ $? -eq 1 ]; then
    if [ $found -eq 0 ]; then
      echo -e "Some configurations were not tested :" 
      found=1
    fi   
    echo -e "   - $testconf not tested for OPENMP"    
  fi  
  mytest=$(ls $testconf/jobcomp_MPI* > /dev/null 2>&1)
  if [ $? -eq 1 ]; then
    if [ $found -eq 0 ]; then
      echo -e "Some configurations were not tested :" 
      found=1
    fi   
    echo -e "   - $testconf not tested for MPI"
  fi  
done

i=1
ierr=0
for testconf in `ls -1 ./Configure_Test/ `;do
  if [ $i -eq 1 ]; then	
  echo "  "
  if [ ${FANCY_OUTPUT} -eq 1 ] ;then
    printf "%26s %14s %20s" COMPILATION EXECUTION REPRODUCIBILITY
  else
    printf "%35s %20s %20s" COMPILATION EXECUTION REPRODUCIBILITY
  fi  
  printf "\n"

  fi
  
  COMPIL_OUT=$(sed '1,1!d' ${testconf}/${testconf}_steps)
  EXEC_OUT=$(sed '2,2!d' ${testconf}/${testconf}_steps)
  REPRO_OUT=$(sed '3,3!d' ${testconf}/${testconf}_steps)
  for TEST in "COMPIL_OUT" "EXEC_OUT" "REPRO_OUT"
  do
    key="${TEST}" 
    eval var='$'$key
    if [ "${var}" == 'Y' ]; then
  	  varname="$(echo -e ${TEST}_PR)"
  	  varvalue="$(echo -e ${FMT_GREENBLD}${var}${FMT_ORD})"
  	  eval "$varname=\$varvalue"
    elif 	[ "${var}" == 'N' ]; then
  	  varname="$(echo -e ${TEST}_PR)"
  	  varvalue="$(echo -e ${FMT_REDBLD}${var}${FMT_ORD})"
  	  eval "$varname=\$varvalue"
  	  ierr=$(($ierr+1))
    else
  	  varname="$(echo -e ${TEST}_PR)"
      varvalue="$(echo -e ${FMT_RVERT}${var}${FMT_ORD})"
  	  eval "$varname=\$varvalue"
    fi
  done
  if [ ${FANCY_OUTPUT} -eq 1 ] ;then
    format="%-12s %20s %29s %29s \n"
  else
    format="%-12s %20s %20s %20s \n"
  fi
  printf "$format" $testconf $COMPIL_OUT_PR $EXEC_OUT_PR $REPRO_OUT_PR

  i=$(($i+1)) 
done  


#
echo "  "
if [ $ierr -eq  0 ]; then  
  echo -e ${FMT_GREEN}""
  cat fancy_sucess.txt
  echo -e ""${FMT_ORD}
  exit 0
elif  [ $ierr -le  3 ]; then
  echo -e ${FMT_ORANGE}""
  cat  fancy_almost.txt
  echo -e ""${FMT_ORD}
  exit 1
elif  [ $ierr -le  9 ]; then
  echo -e ${FMT_RED}""
  cat  fancy_failure.txt
  echo -e ""${FMT_ORD}
  exit 1
else
  echo -e ${FMT_RED2}""
  cat  fancy_critical.txt
  echo -e ""${FMT_ORD}
  exit 1
fi



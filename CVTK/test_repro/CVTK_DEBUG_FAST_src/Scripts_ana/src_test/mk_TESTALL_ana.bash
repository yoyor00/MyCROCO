#!/bin/bash

set -e
set -u
#set -x

./git_process.bash

for testconf in `ls -1 ./Configure_Test/ `;do
#for testconf in BASIN SWASH; do
  [ -d $testconf ] && rm -rf 	$testconf 
  echo $(tput setaf 14 ; tput bold)"=============================="$(tput sgr0)
  echo $(tput setaf 14 ; tput bold)"TESTING $testconf :"$(tput sgr0)
  echo $(tput setaf 14 ; tput bold)"=============================="$(tput sgr0)
  rm -rf $testconf
  ./mk_TestDIR_ana.bash $testconf
  echo "  "  
done

i=1
ierr=0
#for testconf in  SWASH BASIN; do  
for testconf in `ls -1 ./Configure_Test/ `;do
  if [ $i -eq 1 ]; then	
  echo "  "
  printf "%22s %12s %12s" COMPILATION EXECUTION REPRODUCIBILITY
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
  	  varvalue="$(echo -e $(tput setaf 2 ; tput bold )${var}$( tput sgr0))"
  	  eval "$varname=\$varvalue"
    elif 	[ "${var}" == 'N' ]; then
  	  varname="$(echo -e ${TEST}_PR)"
  	  varvalue="$(echo -e $(tput setaf 1 ; tput bold  )${var}$( tput sgr0))"
  	  eval "$varname=\$varvalue"
  	  ierr=$(($ierr+1))
    else
  	  varname="$(echo -e ${TEST}_PR)"
      varvalue="$(echo -e $(  tput setab 7 ; tput bold ;   )${var}$( tput sgr0))"
  	  eval "$varname=\$varvalue"
    fi
  done
  format="%-10s %20s %29s %25s \n"
  printf "$format" $testconf $COMPIL_OUT_PR $EXEC_OUT_PR $REPRO_OUT_PR
   
  i=$(($i+1)) 
done  


env
echo
tty

#
echo "  "
if [ $ierr -eq  0 ]; then  
  tput setaf 2
  cat fancy_sucess.txt
  tput sgr0
  exit 0
elif  [ $ierr -le  3 ]; then
  tput setaf 172
  cat  fancy_almost.txt
  tput sgr0
  exit 1
elif  [ $ierr -le  9 ]; then
  tput setaf 1
  cat  fancy_failure.txt
  tput sgr0
  exit 1
else
  tput setaf 9
  cat  fancy_critical.txt
  tput sgr0
  exit 1
fi



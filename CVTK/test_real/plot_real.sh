#!/bin/bash
set -x

LIST_EXAMPLE=$1
ROOTDIR=$2
NUMBER=$3

if [ ${#ROOTDIR} -eq 0 ]; then
  ROOTDIR=$(dirname $(dirname $PWD))
fi 
if [ ${#LIST_EXAMPLE} -eq 0 ]; then
  LIST_EXAMPLE='BASIN CANYON_A CANYON_B EQUATOR GRAV_ADJ INNERSHELF OVERFLOW SEAMOUNT SHELFRONT SOLITON UPWELLING VORTEX JET RIP  SHOREFACE SWASH THACKER TANK'
fi  

[ ! -d TEST_CASES ] && \cp -r ${ROOTDIR}/Run/TEST_CASES .

REQUIRE="matlab pdfcrop gs"
for i in $REQUIRE
do 
 has_it=$(which $i)
 if [ -z $has_it ]; then
 echo -e "\033[1;31m $i NOT available ... We quit \033[0m" && exit 1
 fi
done

i=0
for EXAMPLE in $LIST_EXAMPLE
  do 
    i=${NUMBER:-$((i=$i+1))}
    echo "-------"
    echo $i
    echo "-------"
    example=$EXAMPLE
    [ "$EXAMPLE" != "IGW" ] && example=$(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')
    myscript="plot_${example}"
    sed -i .bak "s/makepdf\(.*\)=\(.*\)0\(.*\)/makepdf=1/g" TEST_CASES/${myscript}.m
     
    \rm $(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')*.pdf
    \rm $(echo $EXAMPLE |tr '[:lower:]' '[:upper:]')*.pdf
    FILE1=$( ls *.pdf )
    matlab -nodesktop  -nosplash -nodisplay -r "addpath ./TEST_CASES; ${myscript};exit"
    FILE2=$( comm -3 <( ls *.pdf ) <( echo "$FILE1" ) )

    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFIXEDMEDIA -sPAPERSIZE=a4 -dPSFitPage -sOutputFile=tmp.pdf	${FILE2} 

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFIXEDMEDIA -sPAPERSIZE=a4 -dPSFitPage -dSubsetFonts=true -dEmbedAllFonts=true -dPDFSETTINGS=/default  \
-sOutputFile=${EXAMPLE}.pdf                          \
-c "<< \
/EndPage   \
   {    \
      2 eq { pop false }    \
      {    gsave   \
            /Times-Roman 10 selectfont \
         0 -30 moveto ( $( cd $ROOTDIR && git describe --all --long) ) show  \
     /Times-Roman 10 selectfont \
         470 -30 moveto ( $(date "+DATE: %d-%m-%Y TIME: %H:%M:%S") ) show \
         /Times-Bold 18 selectfont  \
         -15 860 moveto ( ${EXAMPLE} ) show  \
         grestore    \
   true   \
} ifelse \
}  bind   \
>> setpagedevice"  -c "<<  /BeginPage  {  0.9 0.9 scale 29.75 42.1 translate  } >> setpagedevice"  -f tmp.pdf 

FILE3=${EXAMPLE}.pdf 
[ $i -eq 1 ] && \rm merged.pdf
[ $i -ne 1 ] && FILE3="merged.pdf $FILE3"
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dSubsetFonts=true -dEmbedAllFonts=true -dPDFSETTINGS=/default  \
-sOutputFile=merged_tmp.pdf ${FILE3}  
\mv merged_tmp.pdf merged.pdf
\rm ${EXAMPLE}.pdf ${FILE2} tmp.pdf

  done


!#/bin/bash

cat QI_LKMODEL.dat | sed 's/(/ /g' | sed 's/)/ /g' | sed 's/,/+i/g' | sed 's/+i-/-i/g' | tr -s ' '| cut -d ' ' -f 1,"$1"| sed 's/ /\t/g'

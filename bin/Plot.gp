#!/usr/bin/env gnuplot
system 'mkdir -p out'
set output 'out/BAND_LKMODEL.pdf'
set terminal pdf
set xtics auto
load "LABELS_LKMODEL.dat"
set ylabel "Energy (eV)"
plot 'BAND_LKMODEL.dat' w line title "LKMODEL band structure"

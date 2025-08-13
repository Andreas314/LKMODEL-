#!/usr/bin/env gnuplot
system 'mkdir -p out'
set output 'out/BAND_LKMODEL.pdf'
set terminal pdf
set xtics auto
load "out/LABELS_LKMODEL.dat"
set ylabel "Energy (eV)"
plot 'out/BAND_LKMODEL.dat' w line title "LKMODEL band structure"

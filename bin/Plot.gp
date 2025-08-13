#!/usr/bin/env gnuplot

set output 'BAND_LKMODEL.pdf'
set terminal pdf
set xtics auto
load "LABELS_LKMODEL.dat"
set ylabel "Energy (eV)"
plot 'BAND_LKMODEL.dat' w line title "LKMODEL band structure"

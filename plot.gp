reset
set grid
set xrange [-10:170]
set yrange [-1:13]
set key left
set xlabel 'E [keV]'
set ylabel 'range [mm]'

plot 'pres.txt' u 1:2 t 'Measure', 'pres.txt' u 3:4 t 'Pres. Geant', 'pres.txt' u 5:6 t 'SRIM', 'ranges.txt' using 1:4:5 with errorbars title 'My Geant'
reset
set grid
set xrange [-10:170]
set yrange [-1:13]
set key left
set xlabel 'E [keV]'
set ylabel 'range [mm]'

plot 'pres.txt' u 1:2 t 'Measure', 'pres.txt' u 3:4 t 'Pres. Geant', 'pres.txt' u 5:6 t 'SRIM', 'ranges.txt' u 1:2 t 'My Geant proj', 'ranges.txt' u 1:3 t 'My Geant'
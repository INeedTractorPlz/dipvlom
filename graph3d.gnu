#! /usr/bin/gnuplot -persist


set xrange [-2:2]
set yrange [-2:2]
set zrange [-0.0025:0.003]

splot "NBodyRK_0.dat" w l,"NBodyRK_1.dat" w l,"NBodyRK_2.dat" w l
pause -1

#! /usr/bin/gnuplot -persist

set terminal png

set output 'graph_r.png'
plot 'NBodyRK_0.dat' using 1:14 w l

set output 'graph_q.png'
plot 'NBodyRK_0.dat' using 1:13 w l

set output 'graph_p.png'
plot 'NBodyRK_0.dat' using 1:12 w l

set output 'graph_psi.png'
plot 'NBodyRK_0.dat' using 1:9 w l

set output 'graph_phi.png'
plot 'NBodyRK_0.dat' using 1:10 w l

set output 'graph_tetta.png'
plot 'NBodyRK_0.dat' using 1:11 w l

set output 'graph_distance.png'
plot 'NBodyRK_0.dat' using 1:2 w l



#!/usr/bin/gnuplot

set terminal pngcairo size 700,524 enhanced font 'Verdana,10'

unset key

set xrange [-3:3]
set yrange [-3:3]
set zrange [-0.004:0.004]

system('rm -rf png3d')
system('mkdir -p  png3d')

do for [ii=1:1000:5] {

set output sprintf('png3d/moon%03.0f.png',ii)

#splot 'fort.1' every ::1::ii w l ls 1 lt rgb 'green', \
#'fort.1' every ::ii::ii w p ls 1 pt 7 ps 0.5 lt rgb 'green' , \
#'fort.2' every ::1::ii w l ls 1 lt rgb 'blue', \
#'fort.2' every ::ii::ii w p ls 1 pt 7 ps 0.5 lt rgb 'blue', \
#'fort.3' every ::1::ii w l ls 1 lt rgb 'red', \
#'fort.3' every ::ii::ii w p ls 1 pt 7 ps 0.5 lt rgb 'red'

splot 'NBodyRK_2.dat' every ::ii::ii  using 2:3:4  ls 1 pt 7 ps 2 lt rgb 'yellow', \
'NBodyRK_1.dat' every ::ii::ii   using 2:3:4  ls 1 pt 7 ps 1 lt rgb 'blue', \
'NBodyRK_0.dat' every ::ii::ii   using 3:4:5  ls 1 pt 7 ps 0.5 lt rgb 'gray'
}

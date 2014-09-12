#! gnuplot
reset

set term postscript eps enhanced color solid size 8.0cm, 6.0cm
set outp 'Frozen.eps'

set key;
set grid front
unset grid
set style data l
set multiplot
set logscale x
unset ylabel
set yrange [0:60]

set lmargin at screen 0.4
set rmargin at screen 0.95

unset ytics
set format x '10^%T'
set xtics 10,10,100
set xlabel 'K'
set label 2 'unstable' at 10,20 front;
plot [1:200] 'InfiniteLayer_N60.dat' using ($2==0.0001 & $4>1 ? $4 : 1./0):3 t 'Frozen coefficient' w l lt 1, \
           '' u ($2==0.0001 & $4>1 ? $4 : 1./0):3 not w filledcurves x1 fillstyle transparent solid 0.4 lt 1, \
           'thresholdtc.dat' t 'Spectral radius' w p pt 5 lt 3 ps 2;

set yrange [0:60]
set ylabel 't_*'
unset log x
set format x 
unset key;
set bmargin at screen 0.132
set tmargin at screen 0.951
set ytics
set xtics 0,0.5,1
set lmargin at screen 0.15
set rmargin at screen 0.4
set xlabel ''
unset label 2
set label 1 'stable' at 0.2, 10 front;
plot [0:1] 'Kless1.dat' using 1:2 not w  filledcurves x1 fillstyle transparent solid 0.4, 'thresholdtc.dat' w p pt 5 lt 3 ps 2

unset multiplot
set outp;
set term wxt enha
!epstopdf Frozen.eps

#! gnuplot
set view map;
set log xy;
# set contour;
# set cntrparam level discrete (1.1), (2), (10)
# set label "1.1" at 50,100 front;
# set label "2" at 100,300 front;
# set label "10" at 100,700 front;
set format xy "$10^{%T}$"
set bmargin 0.05
set tmargin 0.0
set lmargin 0.05
set rmargin 0.05
unset key
set xlabel '$t_1$' offset 0, 0.5
set ylabel '$t_2$' offset 0.5, 0
set cbrange [1:*]
load '/home/shreyas/usr/local/gnuplot-colorbrewer/sequential/OrRd.plt'
splot 'growth.dat' u 1:2:($1<$2?(exp($4)):1/0) w pm3d;
replot 'neutral.dat' ev :::0::0 u 1:2:(1000) w l lt -1 lw 2;

# set term postscript eps enhance color size 8cm,6cm 12
set term epslatex color solid size 8.6cm,8cm
set outp 'figure.tex'
replot
set outp
set term wxt enhance


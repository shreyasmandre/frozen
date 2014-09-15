# ! gnuplot
set view map;
set log xy;
set contour;
set cntrparam level discrete log10(1.1), log10(2), log10(10)
set format xy "$10^{%T}$"
set bmargin 0.05
set tmargin 0.0
set lmargin 0.05
set rmargin 0.05
unset key
set label "1.1" at 10,40 front;
set xr [2:25];
set yr [2:25];
unset xtics
# unset ytics
unset colorbox
splot 'growth.dat' u 1:2:($1<$2?log10(exp($4)):1/0) w pm3d;
replot 'neutral.dat' ev :::0::0 u 1:2:(1000) w l lt -1 lw 2;

# set term postscript eps enhance color size 8cm,6cm 12
set term epslatex color solid size 4cm,4cm
set outp 'figure.tex'
replot
set outp
set term wxt enhance


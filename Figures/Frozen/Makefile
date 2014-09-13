Frozen.eps: Frozen.dvi
	dvips -E -o Frozen.eps Frozen.dvi
Frozen.dvi: Frozen.tex figure.tex figure.eps
	latex Frozen.tex
Frozen.pdf: Frozen.eps
	epstopdf Frozen.eps
figure.tex: Shreyas3.gnu Kless1.dat thresholdtc.dat InfiniteLayer_N60.dat
	gnuplot Shreyas3.gnu

GNUFILE=Color.gnu
all: Neutral1.eps
%.eps: %.dvi
	dvips -E -o $@ $<
Neutral1.dvi: Neutral1.tex figure.tex figure.eps
	latex Neutral1.tex
%.pdf: %.eps
	epstopdf $<
figure.tex: ${GNUFILE} growth.dat neutral.dat
	gnuplot ${GNUFILE}

GNUFILE=Color.gnu
all: Neutral1.eps
%.eps: %.dvi
	dvips -E -o $@ $<
%.dvi: %.tex figure.tex figure.eps
	latex $<
%.pdf: %.eps
	epstopdf $<
figure.tex: ${GNUFILE} growth.dat neutral.dat
	gnuplot ${GNUFILE}

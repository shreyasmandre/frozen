GNUFILE=Color.gnu
all: Neutral0.eps
%.eps: %.dvi
	dvips -E -o $@ $<
%.dvi: %.tex figure.tex figure.eps
	latex Neutral0.tex
%.pdf: %.eps
	epstopdf $<
figure.tex: ${GNUFILE} growth.dat neutral.dat
	gnuplot ${GNUFILE}

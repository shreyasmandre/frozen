all: beta.pdf
beta.pdf: beta.eps
	epstopdf beta.eps
beta.eps: beta.gnu touch.file
	gnuplot beta.gnu
touch.file: bigloop.m
	octave -q bigloop.m


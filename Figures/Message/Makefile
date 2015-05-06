CC=gcc-4.4
INCLUDE=-I/home/shreyas/usr/include
CFLAGS=${INCLUDE} -Wall -O3 -march=native -fomit-frame-pointer -funroll-all-loops
#CFLAGS=${INCLUDE} -Wall -g
LIBS=-L/home/shreyas/usr/lib
#LIBS=-L/home/shreyas/usr/local/lib
LDFLAGS=${LIBS} -larpack -llapack -lblas -latlas -lgfortran -lm
#LDFLAGS=${LIBS} -lfourier -lfftw3_threads -lfftw3 -lpthread -lm
porous: porous.o

clean:
	rm *.o porous

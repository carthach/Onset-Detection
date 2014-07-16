INCLUDES=-I/usr/local/include
LFLAGS=-L/usr/local/lib
LIBS=-lsndfile -lfftw3

onsets: spectral_flux.c
	gcc -o spectral_flux spectral_flux.c $(INCLUDES) $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm spectral_flux

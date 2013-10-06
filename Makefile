all: chimax dft-test
chimax: chimax.cc characters.h slint.h
	g++ -O3 -o chimax chimax.cc -pthread -lfftw3
dft-test: dft-test.cc characters.h slint.h
	g++ -O3 -o dft-test dft-test.cc -lfftw3

all: chimax dft-test test-values test-primitive

chimax: chimax.cc characters.h slint.h
	g++ -O3 -o chimax chimax.cc -pthread -lfftw3
dft-test: dft-test.cc characters.h slint.h
	g++ -O3 -o dft-test dft-test.cc -lfftw3
test-primitive: characters.h slint.h test-primitive.cc
	g++ -O3 -o test-primitive test-primitive.cc -lfftw3
test-values: characters.h slint.h test-values.cc
	g++ -O3 -o test-values test-values.cc -lfftw3

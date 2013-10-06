all: chimax
chimax: chimax.cc characters.h slint.h
	g++ -O3 -o chimax chimax.cc -pthread

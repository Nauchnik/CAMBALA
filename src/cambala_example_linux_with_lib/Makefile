CPP = g++
CPPFLAGS = -O3 -std=c++17 -static

normal_modes: main.o libnormal_modes.a
	${CPP} ${CPPFLAGS} main.o -L. -lnormal_modes -o example_normal_modes

main.o: main.cpp libnormal_modes.a
	${CPP} ${CPPFLAGS} main.cpp -c

clean:
	rm -rf *.o
	rm normal_modes
	clear
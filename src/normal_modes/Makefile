CPP = g++
CPPFLAGS = -O3 -std=c++17 -static

normal_modes: example.o libnormal_modes.a
	${CPP} ${CPPFLAGS} example.o -L. -lnormal_modes -o example_normal_modes

example.o: example.cpp libnormal_modes.a
	${CPP} ${CPPFLAGS} example.cpp -c

clean:
	rm -rf *.o
	rm normal_modes
	clear
	
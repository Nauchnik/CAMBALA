CPP = g++
INC1 = ../alglib/
INCDIRS = -I${INC1}
CPPFLAGS = -Wall -O0 -std=c++11 ${INCDIRS} -g -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS 

pd-sat: main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o sspemdd_parallel.o
	${CPP} ${CPPFLAGS} main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o sspemdd_parallel.o -o SSPEMDD_parallel

specialfunctions.o: ../alglib/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} ../alglib/specialfunctions.cpp -c

linalg.o: ../alglib/linalg.cpp 
	${CPP} ${CPPFLAGS} ../alglib/linalg.cpp -c

ap.o: ../alglib/ap.cpp 
	${CPP} ${CPPFLAGS} ../alglib/ap.cpp -c

alglibmisc.o: ../alglib/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} ../alglib/alglibmisc.cpp -c

alglibinternal.o: ../alglib/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} ../alglib/alglibinternal.cpp -c

sspemdd_parallel.o: sspemdd_parallel.cpp
	${CPP} ${CPPFLAGS} sspemdd_parallel.cpp -c

sspemdd_sequential.o: sspemdd_sequential.cpp
	${CPP} ${CPPFLAGS} sspemdd_sequential.cpp -c

main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c
	
clean:
	rm -rf *.o
	rm SSPEMDD_parallel
	clear

bisect_cpu.o: bisect_cpu.cpp
	${CPP} ${CPPFLAGS} -Wall bisect_cpu.cpp -c

bisect_test.o: bisect_test.cpp
	${CPP} ${CPPFLAGS} -Wall bisect_test.cpp -c

bisect_test: bisect_test.o bisect_cpu.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o
	${CPP} ${CPPFLAGS} bisect_test.o bisect_cpu.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o -o bisect_test

CPP = mpicxx
ALGLIBDIR = ../../../alglib/
INCDIRS = -I${ALGLIBDIR}
CPPFLAGS = -O3 -std=c++0x ${INCDIRS} -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D _MPI

pd-sat: main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o parallel.o
	${CPP} ${CPPFLAGS} main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o parallel.o -o cambala_parallel

specialfunctions.o: ${ALGLIBDIR}specialfunctions.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}specialfunctions.cpp -c

linalg.o: ${ALGLIBDIR}linalg.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}linalg.cpp -c

ap.o: ${ALGLIBDIR}ap.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}ap.cpp -c

alglibmisc.o: ${ALGLIBDIR}alglibmisc.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}alglibmisc.cpp -c

alglibinternal.o: ${ALGLIBDIR}alglibinternal.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}alglibinternal.cpp -c

parallel.o: parallel.cpp
	${CPP} ${CPPFLAGS} parallel.cpp -c

sequential.o: sequential.cpp
	${CPP} ${CPPFLAGS} sequential.cpp -c

main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c
	
clean:
	rm -rf *.o
	rm cambala_parallel
	clear
CPP = mpicxx
INC1 = ../../alglib/
INCDIRS = -I${INC1}
CPPFLAGS = -O3 -std=c++0x ${INCDIRS} -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D _MPI

pd-sat: main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o parallel.o
	${CPP} ${CPPFLAGS} main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o parallel.o -o cambala_parallel

specialfunctions.o: ../../alglib/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} ../../alglib/specialfunctions.cpp -c

linalg.o: ../../alglib/linalg.cpp 
	${CPP} ${CPPFLAGS} ../../alglib/linalg.cpp -c

ap.o: ../../alglib/ap.cpp 
	${CPP} ${CPPFLAGS} ../../alglib/ap.cpp -c

alglibmisc.o: ../../alglib/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} ../../alglib/alglibmisc.cpp -c

alglibinternal.o: ../../alglib/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} ../../alglib/alglibinternal.cpp -c

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
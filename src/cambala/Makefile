CPP = g++
ALGLIBDIR = ../../../alglib
EIGENVALUESDIR = ../../../Eigenvalues
INCDIRS = -I${ALGLIBDIR} -I${EIGENVALUESDIR}
CPPFLAGS = -O3 -std=c++17 -static ${INCDIRS}

cambala: main.o sequential.o normal_modes.a
	${CPP} ${CPPFLAGS} main.o sequential.o normal_modes.a -o cambala

main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c

sequential.o: sequential.cpp
	${CPP} ${CPPFLAGS} sequential.cpp -c
	
normal_modes.o: normal_modes.cpp
	${CPP} ${CPPFLAGS} normal_modes.cpp -c
	
specialfunctions.o: ${ALGLIBDIR}/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}/specialfunctions.cpp -c

linalg.o: ${ALGLIBDIR}/linalg.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}/linalg.cpp -c

ap.o: ${ALGLIBDIR}/ap.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}/ap.cpp -c

alglibmisc.o: ${ALGLIBDIR}/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}/alglibmisc.cpp -c

alglibinternal.o: ${ALGLIBDIR}/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} ${ALGLIBDIR}/alglibinternal.cpp -c

normal_modes.a: alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o normal_modes.o
	ar rcs normal_modes.a alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o normal_modes.o

clean:
	rm -rf *.o
	rm cambala
	rm normal_modes.a
	clear
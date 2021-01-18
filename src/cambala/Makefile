CPP = g++
ALGLIBDIR = ../../../alglib
CPPFLAGS = -O3 -std=c++17 -static

cambala: main.o sequential.o
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} main.o sequential.o -o cambala

main.o: main.cpp
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} main.cpp -c

sequential.o: sequential.cpp
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} sequential.cpp -c
	
specialfunctions.o: ${ALGLIBDIR}/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} ${ALGLIBDIR}/specialfunctions.cpp -c

linalg.o: ${ALGLIBDIR}/linalg.cpp 
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} ${ALGLIBDIR}/linalg.cpp -c

ap.o: ${ALGLIBDIR}/ap.cpp 
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} ${ALGLIBDIR}/ap.cpp -c

alglibmisc.o: ${ALGLIBDIR}/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} ${ALGLIBDIR}/alglibmisc.cpp -c

alglibinternal.o: ${ALGLIBDIR}/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} -I${ALGLIBDIR} ${ALGLIBDIR}/alglibinternal.cpp -c

clean:
	rm -rf *.o
	rm cambala
	clear
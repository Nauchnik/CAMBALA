CPP = g++
INC1 = ../../alglib/
INC2 = ../cambala/
INCDIRS = -I${INC1} -I${INC2}
CPPFLAGS = -O3 -std=c++0x ${INCDIRS}

cambala_smac: cambala_smac.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o
	${CPP} ${CPPFLAGS} cambala_smac.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o -o cambala_smac

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

sequential.o: ../cambala/sequential.cpp
	${CPP} ${CPPFLAGS} ../cambala/sequential.cpp -c

cambala_smac.o: cambala_smac.cpp
	${CPP} ${CPPFLAGS} cambala_smac.cpp -c
	
clean:
	rm -rf *.o
	rm cambala_smac
	clear
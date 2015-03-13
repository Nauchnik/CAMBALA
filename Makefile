CPP = g++
INC1 = ../alglib/
INCDIRS = -I${INC1}
CPPFLAGS = -O3 -std=c++0x ${INCDIRS}

pd-sat: main.o alglibinternal.o alglibmisc.o ap.o integration.o linalg.o specialfunctions.o
	${CPP} ${CPPFLAGS} main.o alglibinternal.o alglibmisc.o ap.o integration.o linalg.o specialfunctions.o -o acoustics_client_app
	
specialfunctions.o: ../alglib/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} ../alglib/specialfunctions.cpp -c

linalg.o: ../alglib/linalg.cpp 
	${CPP} ${CPPFLAGS} ../alglib/linalg.cpp -c

integration.o: ../alglib/integration.cpp 
	${CPP} ${CPPFLAGS} ../alglib/integration.cpp -c

ap.o: ../alglib/ap.cpp 
	${CPP} ${CPPFLAGS} ../alglib/ap.cpp -c

alglibmisc.o: ../alglib/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} ../alglib/alglibmisc.cpp -c

alglibinternal.o: ../alglib/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} ../alglib/alglibinternal.cpp -c

main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c
	
clean:
	rm -rf *.o
	rm acoustics_client_app
	clear
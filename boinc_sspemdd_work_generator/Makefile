SSPEMDD = ..
ALGLIB = ../../alglib

CPP = g++
CPPFLAGS = -O3 \
    -std=c++0x -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
    -DAPP_GRAPHICS -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D NDEBUG \
    -I$(SSPEMDD) \
    -I$(ALGLIB) \
    -L /usr/X11R6/lib \
    -L.

release: sspemdd_work_generator

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o sspemdd_work_generator

sspemdd_work_generator: sspemdd_work_generator.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o
	$(CPP) $(CPPFLAGS) sspemdd_work_generator.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o -o sspemdd_work_generator
	
specialfunctions.o: $(ALGLIB)/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/specialfunctions.cpp -c

linalg.o: $(ALGLIB)/linalg.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/linalg.cpp -c

ap.o: $(ALGLIB)/ap.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/ap.cpp -c

alglibmisc.o: $(ALGLIB)/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/alglibmisc.cpp -c

alglibinternal.o: $(ALGLIB)/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/alglibinternal.cpp -c

sspemdd_sequential.o: $(SSPEMDD)/sspemdd_sequential.cpp
	${CPP} ${CPPFLAGS} $(SSPEMDD)/sspemdd_sequential.cpp -c

sspemdd_work_generator.o: sspemdd_work_generator.cpp
	${CPP} ${CPPFLAGS} sspemdd_work_generator.cpp -c

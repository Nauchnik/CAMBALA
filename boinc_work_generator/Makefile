CAMBALA = ../cambala/
ALGLIB = ../../alglib
MYSQL = /usr/include/mysql
UUID = /usr/include/uuid
INCDIRS  = -I${MYSQL} -I${UUID} -I${ALGLIB} -I${CAMBALA}

CPP_DEBUG_FLAGS = -O0 -std=c++11 -DDEBUG_ALL -g ${INCDIRS}
CPP_RELEASE_FLAGS = -O2 -std=c++11 -D NDEBUG -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS ${INCDIRS} -D _BOINC
MYSQL_LIBS = -luuid /usr/lib/x86_64-linux-gnu/libmysqlclient.so

CPP = g++
CPPFLAGS = $(CPP_RELEASE_FLAGS)

release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: cambala_work_generator

debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
debug: cambala_work_generator

cambala_work_generator: cambala_work_generator.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o
	$(CPP) $(CPPFLAGS) cambala_work_generator.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sequential.o -o cambala_work_generator $(MYSQL_LIBS)
	
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

sequential.o: $(CAMBALA)/sequential.cpp
	${CPP} ${CPPFLAGS} $(CAMBALA)/sequential.cpp -c

cambala_work_generator.o: cambala_work_generator.cpp
	${CPP} ${CPPFLAGS} cambala_work_generator.cpp -c
	
clean:
	rm -fr *.o
	rm cambala_work_generator
	clear

SSPEMDD = ../
ALGLIB = ../../alglib/
BOINC_DIR = ../../../
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2

CPP = g++
CPPFLAGS = -O3 \
    -std=c++0x -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
    -DAPP_GRAPHICS -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D NDEBUG \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -I$(BOINC_ZIP_DIR) \
    -I$(FREETYPE_DIR) \
    -I$(SSPEMDD) \
	-I$(ALGLIB) \
    -L /usr/X11R6/lib \
    -L.

release: sspemdd_boinc

libstdc++.a:
	ln -s `g++ -print-file-name=libstdc++.a`

ttfont.cpp:
	ln -s ../../api/ttfont.cpp .

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o libstdc++.a

sspemdd_boinc: sspemdd_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a
	$(CPP) $(CPPFLAGS) sspemdd_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o -o sspemdd_boinc libstdc++.a -pthread \
	$(BOINC_API_DIR)/libboinc_api.a \
	$(BOINC_LIB_DIR)/libboinc.a

sspemdd_boinc: sspemdd_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o
	${CPP} ${CPPFLAGS} sspemdd_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o -o SSPEMDD
	
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

sspemdd_sequential.o: sspemdd_sequential.cpp
	${CPP} ${CPPFLAGS} sspemdd_sequential.cpp -c

sspemdd_boinc.o: sspemdd_boinc.cpp
	${CPP} ${CPPFLAGS} sspemdd_boinc.cpp -c

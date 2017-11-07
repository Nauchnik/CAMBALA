TARGET := libcambala.a

SOURCES := \
	cambala.cpp \
	scenario.cpp\
	utils.cpp\
	rng.cpp\
	residual/interface.cpp\
	solvers/interface.cpp\
	solvers/discrete.cpp\
	solvers/hillclimbing.cpp\
	solvers/bruteforce.cpp

CXXFLAGS := -g -std=c++11 -O3 

SRC_INCDIRS := ./ 

TGT_LDLIBS  := -leasylogging++
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libeasylogging++.a



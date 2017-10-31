TARGET := libcambala.a

SOURCES := \
	cambala.cpp \
	scenario.cpp\
	utils.cpp\
	residual/selector.cpp\
	residual/cpu64.cpp\
	solvers/discrete.cpp\
	solvers/bruteforce.cpp

CXXFLAGS := -g -std=c++11 -O3 

SRC_INCDIRS := ./ 

TGT_LDLIBS  := -leasylogging++
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libeasylogging++.a



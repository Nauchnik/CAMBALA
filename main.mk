
ADD_TGT :=
ifeq ($(GPU),enable)
	ADD_TGT=src/libcambala/cuda.mk
endif



CXXFLAGS := -g -O3 -Wall -pipe -std=c++11

BUILD_DIR := build
TARGET_DIR := out
INCDIRS = ext/easyloggingpp ext/bbox/COMPI/ ext/bbox/snowgoose/  ext/easyloggingpp/  ext/bbox/LOCSEARCH/ ext/bbox/BBSEARCH/

SUBMAKEFILES := ext/easyloggingpp/easylogging++.mk src/libcambala/libcambala.mk src/cambala/cambala.mk $(ADD_TGT) #src/tests/rescalc.mk

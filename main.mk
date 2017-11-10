CXXFLAGS := -g -O3 -Wall -pipe -std=c++11

BUILD_DIR := build
TARGET_DIR := out
#INCDIRS = ext/easyloggingpp
INCDIRS = ext/easyloggingpp ext/bbox/COMPI/ ext/bbox/snowgoose/  ext/easyloggingpp/  src/tests_ahw/  ext/bbox/LOCSEARCH/ ext/bbox/BBSEARCH/

SUBMAKEFILES := ext/easyloggingpp/easylogging++.mk src/libcambala/libcambala.mk src/cambala/cambala.mk #src/tests/rescalc.mk

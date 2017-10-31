CXXFLAGS := -g -O3 -Wall -pipe -std=c++11

BUILD_DIR := build
TARGET_DIR := out
INCDIRS = ext/easyloggingpp
SUBMAKEFILES := ext/easyloggingpp/easylogging++.mk src/libcambala/libcambala.mk src/cambala/cambala.mk

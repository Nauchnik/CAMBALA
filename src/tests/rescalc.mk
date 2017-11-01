TARGET := test_rescalc

SOURCES := \
	rescalc.cpp

TGT_LDLIBS  := -lcambala -leasylogging++
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libcambala.a libeasylogging++.a

SRC_INCDIRS := ./ ../libcambala/


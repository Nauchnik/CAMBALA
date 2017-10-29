TARGET := cambala

SOURCES := \
	main.cpp

TGT_LDLIBS  := -lcambala
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libcambala.a

SRC_INCDIRS := ./ ../libcambala/


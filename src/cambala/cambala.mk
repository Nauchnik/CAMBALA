TARGET := cambala

SOURCES := \
	main.cpp

SRC_CXXFLAGS := -pthread
TGT_LDLIBS  := -lcambala -leasylogging++ -lgpures -lm -lcudart
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libcambala.a libeasylogging++.a libgpures.a

SRC_INCDIRS := ./ ../libcambala/


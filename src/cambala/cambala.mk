ADD_LDLIBS :=
ADD_PREREQS :=
ADD_CXXFLAGS :=
CUDA_LDLIBS := -lgpures -lm -lcudart
CUDA_PREREQS := libgpures.a
CUDA_CXXFLAGS := -DGPU_ENABLE -pthread
ifeq ($(GPU),enable)
	ADD_LDLIBS=$(CUDA_LDLIBS)
	ADD_PREREQS=$(CUDA_PREREQS)
	ADD_CXXFLAGS=$(CUDA_CXXFLAGS)
endif



TARGET := cambala

SOURCES := \
	main.cpp

SRC_CXXFLAGS := $(ADD_CXXFLAGS)
TGT_LDLIBS  := -lcambala -leasylogging++ $(ADD_LDLIBS)
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libcambala.a libeasylogging++.a $(ADD_PREREQS)

SRC_INCDIRS := ./ ../libcambala/


TARGET := libgpures.a

SOURCES := \
	residual/gpu32.cu

NVCC=nvcc
NVCCFLAGS := -g -std=c++11 -ccbin=/usr/bin/g++-4.8 --compiler-options="-fPIC"  -gencode arch=compute_60,code=sm_60  -O3 -D_FORCE_INLINES

SRC_INCDIRS := ./  ../../ext/cuda_sdk/

TGT_LDLIBS  := -leasylogging++
TGT_LDFLAGS := -L${TARGET_DIR}
TGT_PREREQS := libeasylogging++.a

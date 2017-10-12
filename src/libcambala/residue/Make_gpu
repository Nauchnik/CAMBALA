CPP = g++-4.8
CC = gcc-4.8
CXX = g++-4.8
NVCC = nvcc
INC1 = ../alglib/
INCDIRS = -I${INC1}
CPPFLAGS =  -O3 -Wall -std=c++11 ${INCDIRS} -g -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -fPIC -DFTYPE=${FTYPE}
NVCCFLAGS =  -ccbin=/usr/bin/g++-4.8 --compiler-options="-fPIC"  -gencode arch=compute_60,code=sm_60  -O3 -D_FORCE_INLINES -g -std=c++11 --ptxas-options="-v -dlcm=cg" -I./ ${INCDIRS} -DFTYPE=${FTYPE}

ils_gpu.o: ils_gpu.cu
	${NVCC} ${NVCCFLAGS} ils_gpu.cu -c

ils_cpu.o: ils_cpu.cpp
	${CXX} ${CPPFLAGS} ils_cpu.cpp -c

single: FTYPE = float
single: ils_gpu.o ils_cpu.o

double: FTYPE = double 
double: ils_gpu.o ils_cpu.o

clean:
	rm -rf *.o
	clear

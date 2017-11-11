ROOT = ..
BINS = acoustics_bottom_R_weighted.exe acoustics_bottom_R_uniform.exe acoustics_hydro_R_weighted.exe  acoustics_hydro_R_uniform.exe test_bottom_R_weighted.exe test_bottom_R_uniform.exe test_hydro_R_weighted.exe  test_hydro_R_uniform.exe
TESTS = test_bottom_R_weighted.exe test_bottom_R_uniform.exe test_hydro_R_weighted.exe test_hydro_R_uniform.exe
#LIBS =   $(ROOT)/sspemdd/sspemdd.a $(ROOT)/alglib/alglib.a  $(ROOT)/gpu/ils_cpu.o $(ROOT)/gpu/ils_gpu.o 
LIBS =   $(ROOT)/sspemdd/sspemdd.a $(ROOT)/alglib/alglib.a  
ADDIT_INC = -I$(ROOT)/alglib -I$(ROOT)/sspemdd 

include $(ROOT)/all.inc
-include deps.inc


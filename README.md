SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data

Authors:
- Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS)
- Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
- Vadim Bulavintsev (Matrosov Institute for System Dynamics and Control Theory)



To compile the project execute:
>Make


To compile with CUDA support:
>Make GPU=enable
(on some systems it may be required to edit src/libcambala/cuda.mk and
 uncomment "-ccbin=your_g++_version" to direct to NVCC-compatible G++
 version)

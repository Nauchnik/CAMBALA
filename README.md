CAMBALA: Coupled Acoustic Modes with Bottom Attenuation in Linear Acoustics

Authors:
Pavel Petrov, Andrey Tyshchenko (Il'ichev Pacific Oceanological Institute FEB RAS) 
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory SB RAS)


================

DESCRIPTION

CAMBALA is aimed at solving computationally hard inverse problems in underwater
acoustics. In particular, it can reconstruct the sound speed profile or seabed
parameters in a shallow-water waveguide using a dispersion-based geoacoustic 
inversion scheme.


================

HOW TO BUILD

> Get sources of CAMBALA https://github.com/Nauchnik/Acoustics-at-home/archive/master.zip

> Suppose that the sources are in the folder \Acoustics-at-home

> Download a zip-archive with the last version of the ALGLIB C++ sources from http://www.alglib.net/download.php

> Unzip the archive to a default folder, e.g. \alglib-3.17.0.cpp.gpl. 

> Make a new folder \alglib on the same level as \Acoustics-at-home, copy to it all files from the
subfolder \alglib-3.17.0.cpp.gpl\cpp\src

> Make a new folder \Eigenvalues on the same level as \Acoustics-at-home

> Download sources of Spectra from https://github.com/yixuan/spectra.git

> Copy the folder \Spectra from Spectra\include to the created folder \Eigenvalues

> Download sources of Eigen from https://gitlab.com/libeigen/eigen.git

> Copy the folder \Eigen from \eigen to the created folder \Eigenvalues

Windows:

To build the normal_modes library, use 

- > \Acoustics-at-home\src\VS_normal_modes_static\VS_normal_modes_static.sln

To build CAMBALA.exe, use 

- > \Acoustics-at-home\src\VS_normal_modes\VS_normal_modes.sln

Linux:

To build the normal_modes library:

- > cd ./Acoustics-at-home/src/normal_modes/

- > make

To build the CAMBALA binary:

- > cd ./Acoustics-at-home/src/cambala/

- > make


================

HOW TO CREATE INVERTION SCENARIOS

A manual on creating inversion scenarios is located in
\doc\Creating_an_inversion_scenario_file.pdf


================

HOW TO RUN

> CAMBALA.exe scenario_file_name verbosity

Possible verbosity values: 0, 1, 2.

NB! All additional files specified in a scenario file must be in the folder where a program is launched. 


================

EXAMPLES

Example 1. Invert bottom parameters and the range correction.

> Copy CAMBALA.exe to the folder \scenarios\Example_1_range_cb_rhob\

> CAMBALA.exe scenario1 1

Example 2. Invert bottom parameters and the sound speed profile.

> Copy CAMBALA.exe to the folder \scenarios\Example_2_cb_rhob_ssp\

> CAMBALA.exe scenario2 1



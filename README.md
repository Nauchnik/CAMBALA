CAMBALA: Coupled Acoustic Modes

Authors:
Pavel Petrov (Il'ichev Pacific Oceanological Institute FEB RAS, Far Eastern Federal University) 
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory SB RAS, ITMO University)
Andrey Tyshchenko (Il'ichev Pacific Oceanological Institute FEB RAS, Far Eastern Federal University)

================
DESCRIPTION

CAMBALA is aimed at solving computationally hard inverse problems in underwater
acoustics. In particular, it can reconstruct the sound speed profile or seabed
parameters in a shallow-water waveguide using a dispersion-based geoacoustic 
inversion scheme.

================
HOW TO BUILD

Windows (only sequential version):

> Download CAMBALA zip-archive https://github.com/Nauchnik/Acoustics-at-home/archive/master.zip

> Unzip the archive to the folder \Acoustics-at-home

> Download a zip-archive with ALGLIB C++ sources from http://www.alglib.net/download.php

> Unzip the archive to a default folder, e.g. \alglib-3.17.0.cpp.gpl. 

> make a new folder \alglib, copy to it all files from the subfolder \alglib-3.17.0.cpp.gpl\cpp\src\

> NB! Folders \Acoustics-at-home and \alglib should be on the same level

> Launch a Visual Studio 2019 solution \Acoustics-at-home\src\VS_cambala\VS_cambala.sln

> Compile and build an .exe file in Visual Studio (in the Release mode)

> VS_cambala.exe will appear either in \Acoustics-at-home\src\VS_cambala\x64\Release\

Linux (only MPI version):

- > git clone https://github.com/Nauchnik/Acoustics-at-home.git

- > cd /src/cambala

- > make

================
HOW TO CREATE INVERTION SCENARIOS

A manual on creating inversion scenarios is located in
\doc\Creating_an_inversion_scenario_file.pdf

================
HOW TO RUN

> VS_cambala.exe scenario_file_name verbosity

Possible verbosity values: 0, 1, 2.

NB! All additional files specified in a scenario file must be in the folder where a program is launched. 

================
EXAMPLES

Example 1. Invert bottom parameters and the range correction.

> Copy VS_cambala.exe to the folder \scenarios\Example_1_range_cb_rhob\

> VS_cambala.exe scenario1 1

Example 2. Invert bottom parameters and the sound speed profile.

> Copy VS_cambala.exe to the folder \scenarios\Example_2_cb_rhob_ssp\

> VS_cambala.exe scenario2 1

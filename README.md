CAMBALA: Coupled Acoustic Modes

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

Windows:

> Download sources of from CAMBALA https://github.com/Nauchnik/Acoustics-at-home/archive/master.zip

> Unzip the sources to the folder \Acoustics-at-home

> Download a zip-archive with ALGLIB C++ sources from http://www.alglib.net/download.php

> Unzip the archive to a default folder, e.g. \alglib-3.17.0.cpp.gpl. 

> Make a new folder \alglib on the same level as \Acoustics-at-home, copy to it all files from the
subfolder \alglib-3.17.0.cpp.gpl\cpp\src

> Make a new folder \Eigenvalues on the same level as \Acoustics-at-home

> Download sources of Spectra from https://github.com/yixuan/spectra.git

> Unzip the sources, copy the folder \include\Spectra to the created folder \Eigenvalues

> Open a Visual Studio solution \Acoustics-at-home\src\VS_cambala\VS_cambala.sln

> Compile and build an .exe file in Visual Studio (in the Release mode)

> CAMBALA.exe will appear, e.g., in \Acoustics-at-home\src\VS_cambala\x64\Release\

Linux:

- > git clone https://github.com/Nauchnik/Acoustics-at-home.git

- > cd /src/cambala/

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

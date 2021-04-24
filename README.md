CAMBALA: Coupled Acoustic Modes with Bottom Attenuation in Linear Acoustics

Authors:
Pavel Petrov*, Oleg Zaikin**, Andrey Tyshchenko*
* Il'ichev Pacific Oceanological Institute FEB RAS
** Matrosov Institute for System Dynamics and Control Theory SB RAS


================

DESCRIPTION

CAMBALA is aimed at computing sound propagation in 3D-varying shallow water waveguides.

CAMBALA Inversion Tools is aimed at solving computationally hard inverse problems
in underwater acoustics. In particular, it can reconstruct the sound speed profile
or seabed parameters in a shallow water waveguide using a dispersion-based geoacoustic 
inversion scheme.


================

HOW TO BUILD

> Get sources of CAMBALA https://github.com/Nauchnik/CAMBALA/archive/master.zip

> Suppose that the sources are in the directory \CAMBALA

> Download a zip-archive with the last version of the ALGLIB C++ sources from http://www.alglib.net/download.php

> Unzip the archive to a default directory, e.g. \alglib-3.17.0.cpp.gpl. 

> Make a new directory \alglib on the same level as \CAMBALA, copy to it all files from the
subdirectory \alglib-3.17.0.cpp.gpl\cpp\src

> Make a new directory \Eigenvalues on the same level as \CAMBALA

> Download sources of Spectra from https://github.com/yixuan/spectra.git

> Copy the directory \Spectra from Spectra\include to the created directory \Eigenvalues

> Download sources of Eigen from https://gitlab.com/libeigen/eigen.git

> Copy the directory \Eigen from \eigen to the created 
\Eigenvalues

Windows:

To build the CAMBALA library, use 

- > \CAMBALA\src\VS_normal_modes_static\VS_normal_modes_static.sln

To build CAMBALA.exe, use 

- > \CAMBALA\src\VS_normal_modes\VS_normal_modes.sln

Linux:

To build the CAMBALA library:

- > cd ./CAMBALA/src/normal_modes/

- > make

To build the CAMBALA binary:

- > cd ./CAMBALA/src/normal_modes/

- > make


================

HOW TO CREATE INVERTION SCENARIOS FOR CAMBALA INVERSION TOOLS

A manual on creating inversion scenarios is located in
\doc\Creating_an_inversion_scenario_file.pdf


================

HOW TO RUN CAMBALA INVERSION TOOLS

> CAMBALA_Inversion_Tools.exe scenario_file_name verbosity

Possible verbosity values: 0, 1, 2.

NB! All additional files specified in a scenario file must be in the directory
where a program is launched. 


================

EXAMPLES for CAMBALA INVERSION TOOLS

Example 1. Invert bottom parameters and the range correction.

> Copy CAMBALA_Inversion_Tools.exe to the directory \scenarios\Example_1_range_cb_rhob\

> CAMBALA_Inversion_Tools.exe scenario1 1

Example 2. Invert bottom parameters and the sound speed profile.

> Copy CAMBALA_Inversion_Tools.exe to the directory \scenarios\Example_2_cb_rhob_ssp\

> CAMBALA_Inversion_Tools.exe scenario2 1



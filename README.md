### CAMBALA: Coupled Acoustic Modes with Bottom Attenuation in Linear Acoustics

Authors:
Pavel Petrov^, Oleg Zaikin^^, Andrey Tyshchenko^

^Il'ichev Pacific Oceanological Institute FEB RAS

^^Matrosov Institute for System Dynamics and Control Theory SB RAS

### DESCRIPTION

CAMBALA is aimed at computing sound propagation in 3D-varying shallow water waveguides.

CAMBALA Inversion tools is aimed at solving computationally hard inverse problems
in underwater acoustics. In particular, it can reconstruct the sound speed profile
or seabed parameters in a shallow water waveguide using a dispersion-based geoacoustic 
inversion scheme.

### HOW TO BUILD

#### CAMBALA

##### Windows

Use the Visual studio solution /src/VS_normal_modes/VS_normal_modes.sln

##### Linux

> git clone --recurse-submodules https://github.com/Nauchnik/CAMBALA.git

> cd src/normal_modes

> cmake -DCMAKE_BUILD_TYPE=Release .

> make

#### CAMBALA Inversion tools

##### Windows

Use the Visual studio solution ./src/VS_inversion_tools/VS_inversion_tools.sln

##### Linux

> git clone --recurse-submodules https://github.com/Nauchnik/CAMBALA.git

> cd src/inversion_tools

> cmake -DCMAKE_BUILD_TYPE=Release .

> make

### HOW TO CREATE INVERTION SCENARIOS FOR CAMBALA INVERSION TOOLS

A manual on creating inversion scenarios is here
\doc\Creating_an_inversion_scenario_file.pdf

### HOW TO RUN CAMBALA INVERSION TOOLS

> CAMBALA_inversion_tools.exe scenario_file_name verbosity

Possible verbosity values: 0, 1, 2.

NB! All additional files specified in a scenario file must be in the directory
where a program is launched. 

### EXAMPLES FOR CAMBALA INVERSION TOOLS

NB! In all examples below CAMBALA_inversion_tools.exe must be first copied to the corresponding directory.

Example 1. Invert bottom parameters and the range correction.
Directory: \scenarios\Example_1_range_cb_rhob\

> CAMBALA_inversion_tools.exe scenario1 1

Example 2. Invert bottom parameters and the sound speed profile.
Directory: \scenarios\Example_2_cb_rhob_ssp\

> CAMBALA_inversion_tools.exe scenario2 1

### CITATION

CAMBALA Inversion tools can be cited as follows:

```
@article{ZaikinP16-OIDP,
author  = "Oleg Zaikin and Pavel Petrov",
title   = "Algorithm of reconstruction of the sound speed profile in a shallow-water geoacoustic waveguide from modal dispersion data",
journal = "Optoelectronics, Instrumentation and Data Processing",
year    = "2016",
volume  = "52",
number  = "3",
pages   = "259--265"
}

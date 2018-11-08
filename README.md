CAMBALA: Coupled Acoustic Modes

Authors:
Pavel Petrov (Il'ichev Pacific Oceanological Institute FEB RAS, Far Eastern Federal University), 
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory SB RAS, ITMO University)

================================================================================
DESCRIPTION

CAMBALA is aimed at solving computationally hard inverse problems in underwater
acoustics. In particular, it can reconstruct the sound speed profile or seabed
parameters in a shallow-water waveguide using a dispersion-based geoacoustic 
inversion scheme.

================================================================================
HOW TO BUILD

Windows (only sequential version):
> Download the zip-archive https://github.com/Nauchnik/Acoustics-at-home/archive/master.zip
> Unzip the archive
> Launch VS_cambala\VS_cambala.sln - it is a Visual Studio 2017 solution
> Compile and build in Visual Studio

Linux (only MPI version):
> git clone https://github.com/Nauchnik/Acoustics-at-home.git
> cd /cambala
> make

================================================================================
HOW TO LAUNCH

> VS_cambala.exe scenario_file_name verbosity

Possible verbosity values: 0, 1, 2.

NB! All additional files specified in a scenario file must be in the folder where a program is launched. 

================================================================================
EXAMPLES

Example 1. Invert bottom parameters and the range correction.

> Copy VS_cambala.exe to the folder \scenarios\Example_1_range_cb_rhob\
> VS_cambala.exe scenario1 1

Example 2. Invert bottom parameters and the sound speed profile.

> Copy VS_cambala.exe to the folder \scenarios\Example_2_cb_rhob_ssp\
> VS_cambala.exe scenario2 1

================================================================================
HOW TO CREATE INVERTION SCENARIOS

A manual on creating inversion scenarios is located in
\doc\Creating_an_inversion_scenario_file.docx



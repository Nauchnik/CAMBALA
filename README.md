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

In a scenario file constant and variable parameters should be specified
Possible verbosity values: 0, 1, 2.
  
NB! all additional files specified in a scenario file must be in the same folder. 

================================================================================
EXAMPLES

> VS_cambala.exe test_scenario.txt 1

test_scenario.txt is located in the \scenarios folder
It contains the following strings:
% constant parameters
dtimes_file 50_ac_modes_R7km_dtimes.txt
spmag_file(string|no) no
function_type uniform
h 50
H 300
cw0 1500
cw1 1500
cw2 1498
cw3 1496
cw4 1493
cb 1700
tau 0
% variable parameters
R 6995:5:7005
rhob 1.6:0.1:1.8

Comment. Here the file 50_ac_modes_R7km_dtimes.txt is specified, so it should be in the same folder.
Two paremeters are variable: R and rhob, left and right bounds as well as the step are specified for them
in the format left_bound:step:right_bound. Values of all other parameters are constanst.
Both R and rhob have 3 possible values, that gives us the search space of 9 points.



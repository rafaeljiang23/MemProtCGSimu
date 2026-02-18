# MemProtCGSimu
This is a workflow to simulate membrane protein (TRPV3) oligomer plasiticity through membrane-diffusive protomer exchange via a patchy particle model and coarse-grained molecular dynamics simulation. The codes are developed in BIO-AFM-LAB at Weill Cornell Medicine under the supervision of Professor Simon Scheuring.

Developer: Yining Jiang

Publication: XXX (NOTE: This repo is currently made public solely for peer-review)

User should email to the corresponding author of the paper for details about this work: Professor Simon Scheuring (sis2019@med.cornell.edu)

NOTE: Any usage of the codes should cite the publication mentioned above.

## System requirements:
1. Operating system for code development : macOS Big Sur Version 11.7.8
2. Software for code development: MATLAB (MathWorks) 2023b
3. Additional add-ons: MIJI
4. Non-standard hardware: N/A

## Installation instructions: 
1. The codes require installation of MATLAB (MathWorks) 2023b. An installation guide can be found at: https://www.mathworks.com/help/install/.
2. MIJI is recommanded (not required) for visualizing data. An installation guide can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab.
3. The installation should take less than one hour on a "normal" desktop computer

## General instructions:
### Source code
#### Main scripts:
1. mem_simu_v10b.m  (core simulation codes)
2. run_MemCGSimu_Assembly.m   (master script)
3. create_EE_transition_geo.m   (codes for constructing the energy landscape matrices)
4. analyze_v5.m  (codes for analyzing the simulation results)
5. geometry_trpv3_v2.m   (codes for defining TRPV3 protomer geometric descriptor)
#### Helper functions:
1. analyze_v3_init.m
2. supplement.m

#### Instruction
User should work on the master script, and change parameters in the core simulation codes with caution. 
A pre-defined simulation setup file (.mat) is required to define the intial (x,y,a) position and number of particles, the simulation field size, bond etc. A test .mat file is provided in the 'test/initial setup' directory (user should download this file to the same directory of the source codes).

### Test
#### initial setup
1. init_coor_geo_v2_N64_L32.mat (pre-defined setup file for simulation with ~50% protein coverage)

#### display results
1. display_transition_tetra_N64_L32_geo_v2_test63c_1.mat
2. display2_transition_tetra_N64_L32_geo_v2_test63c_1.mat

#### Instruction
For setting up the simulation, user should download 'init_coor_geo_v2_N64_L32.mat' to the same directory of the source codes
For displaying the simulation results, user should refer to the master script.

## Demo (Test data)
A test data is provided in 'test' directory. Two outcome files are provided ('display_transition_tetra_N64_L32_geo_v2_test63c_1.mat' and 'display2_transition_tetra_N64_L32_geo_v2_test63c_1.mat') to give main outcome files from the test simulation (100s, setup file: 'init_coor_geo_v2_N64_L32.mat'). Note that the simulation codes provided here is designed for coarse-grained molecular dynamics simulations with no additional long-range force field between particles, and no interrupting particles (these features are simple add-ons to the codes provided here). The energy landscape matrices were designed for the TRPV3 tetramer/pentamer transitions, with parameter optimized to reflect the equlibrium energy difference between these two stable TRPV3 species.
These tests are expected to run for less than 30 minutes for demo on a "normal" desktop computer following the instructions provided in the main scripts.

## Instruction for use
User should create simulation setup file (.mat), which should be placed in the same directory with other source codes (template: 'init_coor_geo_v2_N64_L32.mat').
User can potentially use their own energy landscape or load pre-defined energy landscape file (.mat) (template: 'init_EE_geo_v2_transition_test63c.mat', created during test run).
User should adjust the parameters in the main scripts accourding to the nature of their simulations. 

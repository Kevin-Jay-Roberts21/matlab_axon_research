% Axon voltage plotter. This code is meant to plot data from the following
% files

% Additionally, variables like list_of_positions, list_of_times, and pause 
% may have to be modified depending on HH, SC or DC to obtain useful animations
% and/or plots.
% 
% Kevin Roberts
% January 2025

% Plotting HH or SC or DC code can be

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%
% DATA USED IN PAPER %
%%%%%%%%%%%%%%%%%%%%%%

% HH_1 = load('projects/axon_simulations/Paper_HH_Data/HH_1.mat');
% HH_2 = load('projects/axon_simulations/Paper_HH_Data/HH_2.mat');
% HH_3 = load('projects/axon_simulations/Paper_HH_Data/HH_3.mat');

% SC_Cohen_set1 = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1.mat');
% SC_Cohen_set1_long = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1_long.mat');
% SC_Cohen_set2 = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set2.mat');
% SC_Cohen_set1_T33 = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1_T33.mat');
% SC_Cohen_set1_T33_long = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1_T33_long.mat');
% SC_Cohen_set1_new_conductances = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1_new_conductances.mat');
% DC_Cohen_set1 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1.mat');
% DC_Cohen_set1_long = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_long.mat');
% DC_Cohen_set2 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set2.mat');
% DC_Cohen_set1_T33 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_T33.mat');
% DC_Cohen_set1_T33_long = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_T33_long.mat');
% DC_Cohen_set1_new_conductances = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_new_conductances.mat');

% SC_Huang_Myelinated_set1 = load('projects/axon_simulations/Paper_Huang_Data/SC_Huang_Myelinated_set1.mat');
% SC_Huang_Tube_set1 = load('projects/axon_simulations/Paper_Huang_Data/SC_Huang_Tube_set1.mat');
% SC_Huang_TubeParalyene_set1 = load('projects/axon_simulations/Paper_Huang_Data/SC_Huang_TubeParalyene_set1.mat');
% DC_Huang_Myelinated_set1 = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Myelinated_set1.mat');
% DC_Huang_Myelinated_set1_rpn10fold = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Myelinated_set1_rpn10fold.mat');
% DC_Huang_Myelinated_set1_rpn100fold = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Myelinated_set1_rpn100fold.mat');
% DC_Huang_Tube_set1 = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Tube_set1.mat');
% DC_Huang_Tube_set1_rpn10fold = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Tube_set1_rpn10fold.mat');
% DC_Huang_Tube_set1_rpn100fold = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_Tube_set1_rpn100fold.mat');
% DC_Huang_TubeParalyene_set1_rpn100fold = load('projects/axon_simulations/Paper_Huang_Data/DC_Huang_TubeParalyene_set1_rpn100fold.mat');

% PARAMETER SENSITIVITY
% DC_Cohen_set1_rpn1 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1.mat');
% DC_Cohen_set1_rpn25 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn25.mat');
% DC_Cohen_set1_rpn50 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn50.mat');
% DC_Cohen_set1_rpn200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn200.mat');
% DC_Cohen_set1_rpn400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn400.mat');
% DC_Cohen_set1_rpn600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn600.mat');
% DC_Cohen_set1_rpn800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn800.mat');
% DC_Cohen_set1_rpn1000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1000.mat');
% DC_Cohen_set1_rpn1200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1200.mat');
% DC_Cohen_set1_rpn1400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1400.mat');
% DC_Cohen_set1_rpn1600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1600.mat');
% DC_Cohen_set1_rpn1800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn1800.mat');
% DC_Cohen_set1_rpn2000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn2000.mat');
% DC_Cohen_set1_rpn2200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn2200.mat');
% DC_Cohen_set1_rpn2400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn2400.mat');
% DC_Cohen_set1_rpn2600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn2600.mat');
% DC_Cohen_set1_rpn2800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn2800.mat');
% DC_Cohen_set1_rpn3000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_rpn3000.mat');
% DC_Cohen_set2_rpn1 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1.mat');
% DC_Cohen_set2_rpn25 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn25.mat');
% DC_Cohen_set2_rpn50 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn50.mat');
% DC_Cohen_set2_rpn200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn200.mat');
% DC_Cohen_set2_rpn400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn400.mat');
% DC_Cohen_set2_rpn600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn600.mat');
% DC_Cohen_set2_rpn800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn800.mat');
% DC_Cohen_set2_rpn1000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1000.mat');
% DC_Cohen_set2_rpn1200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1200.mat');
% DC_Cohen_set2_rpn1400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1400.mat');
% DC_Cohen_set2_rpn1600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1600.mat');
% DC_Cohen_set2_rpn1800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn1800.mat');
% DC_Cohen_set2_rpn2000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn2000.mat');
% DC_Cohen_set2_rpn2200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn2200.mat');
% DC_Cohen_set2_rpn2400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn2400.mat');
% DC_Cohen_set2_rpn2600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn2600.mat');
% DC_Cohen_set2_rpn2800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn2800.mat');
% DC_Cohen_set2_rpn3000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_rpn3000.mat');

% DC_Cohen_set1_Cmy01 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy01.mat');
% DC_Cohen_set1_Cmy02 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy02.mat');
% DC_Cohen_set1_Cmy03 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy03.mat');
% DC_Cohen_set1_Cmy04 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy04.mat');
% DC_Cohen_set1_Cmy05 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy05.mat');
% DC_Cohen_set1_Cmy06 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy06.mat');
% DC_Cohen_set1_Cmy07 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy07.mat');
% DC_Cohen_set1_Cmy08 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy08.mat');
% DC_Cohen_set1_Cmy09 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy09.mat');
% DC_Cohen_set1_Cmy10 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy10.mat');
% DC_Cohen_set1_Cmy11 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy11.mat');
% DC_Cohen_set1_Cmy12 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy12.mat');
% DC_Cohen_set1_Cmy13 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy13.mat');
% DC_Cohen_set1_Cmy14 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy14.mat');
% DC_Cohen_set1_Cmy15 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy15.mat');
% DC_Cohen_set1_Cmy16 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy16.mat');
% DC_Cohen_set1_Cmy17 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy17.mat');
% DC_Cohen_set1_Cmy18 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy18.mat');
% DC_Cohen_set1_Cmy19 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy19.mat');
% DC_Cohen_set1_Cmy20 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy20.mat');
% DC_Cohen_set2_Cmy01 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy01.mat');
% DC_Cohen_set2_Cmy02 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy02.mat');
% DC_Cohen_set2_Cmy03 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy03.mat');
% DC_Cohen_set2_Cmy04 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy04.mat');
% DC_Cohen_set2_Cmy05 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy05.mat');
% DC_Cohen_set2_Cmy06 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy06.mat');
% DC_Cohen_set2_Cmy07 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy07.mat');
% DC_Cohen_set2_Cmy08 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy08.mat');
% DC_Cohen_set2_Cmy09 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy09.mat');
% DC_Cohen_set2_Cmy10 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy10.mat');
% DC_Cohen_set2_Cmy11 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy11.mat');
% DC_Cohen_set2_Cmy12 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy12.mat');
% DC_Cohen_set2_Cmy13 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy13.mat');
% DC_Cohen_set2_Cmy14 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy14.mat');
% DC_Cohen_set2_Cmy15 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy15.mat');
% DC_Cohen_set2_Cmy16 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy16.mat');
% DC_Cohen_set2_Cmy17 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy17.mat');
% DC_Cohen_set2_Cmy18 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy18.mat');
% DC_Cohen_set2_Cmy19 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy19.mat');
% DC_Cohen_set2_Cmy20 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Cmy20.mat');
% 
% DC_Cohen_set1_Rmy50 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy50.mat');
% DC_Cohen_set1_Rmy100 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy100.mat');
% DC_Cohen_set1_Rmy200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy200.mat');
% DC_Cohen_set1_Rmy300 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy300.mat');
% DC_Cohen_set1_Rmy400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy400.mat');
% DC_Cohen_set1_Rmy500 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy500.mat');
% DC_Cohen_set1_Rmy600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy600.mat');
% DC_Cohen_set1_Rmy700 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy700.mat');
% DC_Cohen_set1_Rmy800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy800.mat');
% DC_Cohen_set1_Rmy900 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy900.mat');
% DC_Cohen_set1_Rmy1000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Rmy1000.mat');
% DC_Cohen_set2_Rmy50 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy50.mat');
% DC_Cohen_set2_Rmy100 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy100.mat');
% DC_Cohen_set2_Rmy200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy200.mat');
% DC_Cohen_set2_Rmy300 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy300.mat');
% DC_Cohen_set2_Rmy400 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy400.mat');
% DC_Cohen_set2_Rmy500 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy500.mat');
% DC_Cohen_set2_Rmy600 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy600.mat');
% DC_Cohen_set2_Rmy700 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy700.mat');
% DC_Cohen_set2_Rmy800 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy800.mat');
% DC_Cohen_set2_Rmy900 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy900.mat');
% DC_Cohen_set2_Rmy1000 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Rmy1000.mat');

% DC_Cohen_set1_Ri20 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri20.mat');
% DC_Cohen_set1_Ri40 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri40.mat');
% DC_Cohen_set1_Ri60 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri60.mat');
% DC_Cohen_set1_Ri80 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri80.mat');
% DC_Cohen_set1_Ri100 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri100.mat');
% DC_Cohen_set1_Ri120 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri120.mat');
% DC_Cohen_set1_Ri140 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri140.mat');
% DC_Cohen_set1_Ri160 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri160.mat');
% DC_Cohen_set1_Ri180 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri180.mat');
% DC_Cohen_set1_Ri200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Ri200.mat');
% DC_Cohen_set2_Ri20 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri20.mat');
% DC_Cohen_set2_Ri40 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri40.mat');
% DC_Cohen_set2_Ri60 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri60.mat');
% DC_Cohen_set2_Ri80 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri80.mat');
% DC_Cohen_set2_Ri100 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri100.mat');
% DC_Cohen_set2_Ri120 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri120.mat');
% DC_Cohen_set2_Ri140 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri140.mat');
% DC_Cohen_set2_Ri160 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri160.mat');
% DC_Cohen_set2_Ri180 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri180.mat');
% DC_Cohen_set2_Ri200 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set2_Ri200.mat');

% TEMPERATURE ANALYSIS
% HH_Tbase = load('projects/axon_simulations/Paper_Temp_Data/HH_Tbase.mat');
% HH_T7 = load('projects/axon_simulations/Paper_Temp_Data/HH_T7.mat');
% HH_T9 = load('projects/axon_simulations/Paper_Temp_Data/HH_T9.mat');
% HH_T11 = load('projects/axon_simulations/Paper_Temp_Data/HH_T11.mat');
% HH_T13 = load('projects/axon_simulations/Paper_Temp_Data/HH_T13.mat');
% HH_T15 = load('projects/axon_simulations/Paper_Temp_Data/HH_T15.mat');
% HH_T17 = load('projects/axon_simulations/Paper_Temp_Data/HH_T17.mat');
% HH_T19 = load('projects/axon_simulations/Paper_Temp_Data/HH_T19.mat');
% HH_T21 = load('projects/axon_simulations/Paper_Temp_Data/HH_T21.mat');
% HH_T23 = load('projects/axon_simulations/Paper_Temp_Data/HH_T23.mat');
% HH_T25 = load('projects/axon_simulations/Paper_Temp_Data/HH_T25.mat');
% HH_T27 = load('projects/axon_simulations/Paper_Temp_Data/HH_T27.mat');
% HH_T29 = load('projects/axon_simulations/Paper_Temp_Data/HH_T29.mat');
% HH_T31 = load('projects/axon_simulations/Paper_Temp_Data/HH_T31.mat');
HH_T33 = load('projects/axon_simulations/Paper_Temp_Data/HH_T33.mat'); % decays
HH_T33_diff = load('projects/axon_simulations/Paper_Temp_Data/HH_T33_different_stim_time.mat');
% HH_T35 = load('projects/axon_simulations/Paper_Temp_Data/HH_T35.mat'); % decays
% HH_T37 = load('projects/axon_simulations/Paper_Temp_Data/HH_T37.mat'); % decays

% SC_Cohen_set1_T20 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T20.mat');
% SC_Cohen_set1_T22 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T22.mat');
% SC_Cohen_set1_T24 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T24.mat');
% SC_Cohen_set1_T26 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T26.mat');
% SC_Cohen_set1_T28 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T28.mat');
% SC_Cohen_set1_T30 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T30.mat');
% SC_Cohen_set1_T32 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T32.mat');
% SC_Cohen_set1_T34 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T34.mat');
% SC_Cohen_set1_T36 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T36.mat');
% SC_Cohen_set1_T38 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T38.mat');
% SC_Cohen_set1_T40 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T40.mat');
% SC_Cohen_set1_T42 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T42.mat');
% SC_Cohen_set1_T44 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T44.mat');
% SC_Cohen_set1_T46 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T46.mat');
% SC_Cohen_set1_T48 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T48.mat');
% SC_Cohen_set1_T50 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T50.mat');
% SC_Cohen_set1_T52 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T52.mat');
% SC_Cohen_set1_T54 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T54.mat');
% SC_Cohen_set1_T56 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T56.mat');
SC_Cohen_set1_T58 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T58.mat'); % decays
% SC_Cohen_set1_T60 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T60.mat'); % decays
% SC_Cohen_set1_T62 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T62.mat'); % decays
% SC_Cohen_set1_T64 = load('projects/axon_simulations/Paper_Temp_Data/SC_Cohen_set1_T64.mat'); % decays
% 
% DC_Cohen_set1_T20 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T20.mat');
% DC_Cohen_set1_T22 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T22.mat');
% DC_Cohen_set1_T24 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T24.mat');
% DC_Cohen_set1_T26 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T26.mat');
% DC_Cohen_set1_T28 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T28.mat');
% DC_Cohen_set1_T30 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T30.mat');
% DC_Cohen_set1_T32 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T32.mat');
DC_Cohen_set1_T34 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T34.mat'); % decays
% DC_Cohen_set1_T36 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T36.mat'); % decays
% DC_Cohen_set1_T38 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T38.mat'); % decays
% DC_Cohen_set1_T40 = load('projects/axon_simulations/Paper_Temp_Data/DC_Cohen_set1_T40.mat'); % decays


% picking interval
interval1 = [0.0400 0.0485]; % interval is in cm
interval2 = [0.0800 0.0885]; % interval is in cm
interval3 = [0.1200 0.1285];
time_shot = 2.75; % in ms
space_shot = 0.12; % in cm

% picking time shots
% time1 = 5.1; % in ms
% time2 = 9.5; % in ms
% time3 = 10; % in ms
time1 = 1.1; % in ms
time2 = 3.25; % in ms
time3 = 3.75; % in ms
list_of_times = {time1, time2, time3};

% picking space shots
position1 = 1.1; % in cm
position2 = 2; % in cm
position3 = 3; % in cm
position4 = 4; % in cm
position5 = 5; % in cm

% position1 = 0.0800; % in cm
% position2 = 0.0828; % in cm
% % position3 = 0.0857; % in cm
% position4 = 0.0885; % in cm

% position1 = 0.0005; % in cm
% position2 = 0.04; % in cm
% position3 = 0.08; % in cm
% position4 = 0.12; % in cm
% position5 = 0.16; % in cm

list_of_positions = {position1, position2, position3, position4, position5};


% picking pause (this controls the speed of the animation, the pause variable 
% is in seconds. It is how many seconds each time (or space) shot will be
% paused at). NOTE: the legends are what slows down the animations, may
% condsider getting rid of or adding them in certain cases.
p = 0.001;

% creating a set of data from multiple experiments used to plot animation
% (the first element in the set_of_data is darkred, then the proceeding elements get
% brighter and brighter until the last element which is the brightest red)
% set_of_data = {DC_Cohen_set1_rpn50, DC_Cohen_set1_rpn1, DC_Cohen_set2_rpn50, DC_Cohen_set2_rpn1};
% set_of_data = {DC_Cohen_set1_T30, DC_Cohen_set1_T32, DC_Cohen_set1_T34, DC_Cohen_set1_T36, DC_Cohen_set1_T38, DC_Cohen_set1_T40};

data = HH_T33_diff;
 
% plot_zoomed_in_region_w_AP_at_spaces(data, time_shot, interval1, interval2, interval3);
% plot_Vm_minus_Vmy_picture(data, time_shot);
% plot_animation_voltage_vs_time(data, p);
% plot_animation_voltage_vs_space(data, p);
% plot_animation_probabilities_vs_time(data, p);
% plot_animation_probabilities_vs_space(HH_data_Temp_32, p);
plot_time_and_space_shots(data, list_of_positions, list_of_times);
% plot_voltage_vs_time_comparison_animation(set_of_data, p);
% plot_voltage_vs_space_comparison_animation(set_of_data, p);
% plot_voltage_vs_time_comparison(set_of_data, space_shot);
% plot_voltage_vs_space_comparison(set_of_data, time_shot);
% plot_Vm_and_Vm_minus_Vmy_vs_space(data, p);
% plot_voltage_vs_space_comparison_variable_dt(set_of_data1, p);

%%%%%%%%%%%%%%%%%%%%%
% PLOTTER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%
function plot_zoomed_in_region_w_AP_at_spaces(data, time_shot, interval1, interval2, interval3)
    % Unpack values
    L = data.L;         % Total axon length (cm)
    m = data.m;         % Number of spatial points
    dt = data.dt;       % Time step (ms)
    dx = data.dx;       % Spatial step (cm)

    idx_time = round(time_shot / dt);  % Time index

    % --- Figure 1: Voltage snapshot at a specific time ---
    figure(1)
    t_full = linspace(0, L, m);
    plot(t_full, data.Vm_all(idx_time, :), 'b-');
    hold on
    plot(t_full, data.Vm_minus_Vmy(idx_time, :), 'r-');
    hold off

    axis([0 L -70 30]);
    text(0.125, 12, sprintf('Time: %.3f ms', time_shot), 'FontSize', 9, 'BackgroundColor', 'w');
    legend('$V_m$', '$V_m - V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex');
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel('Length of the axon in cm.')
    title('$V_m$ and $V_m - V_{my}$ for Set (2) on the DC model', 'Interpreter', 'latex');

    % --- Figure 2: Vm vs time at positions in the intervals ---
    figure(2)
    hold on

    % Time vector
    total_timesteps = size(data.Vm_all, 1);
    t_vec = (0:total_timesteps-1) * dt;

    % --- Interval 1: Axon Segment 5 ---
    x_positions1 = linspace(interval1(1), interval1(2), 4);
    spatial_indices1 = round(x_positions1 / dx) + 1;

    % Interval range for legend
    interval1_range = sprintf('%.2f-%.4f cm', interval1(1), interval1(2));

    for i = 1:length(spatial_indices1)
        xi = spatial_indices1(i);
        Vm_trace = data.Vm_all(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'b-', 'DisplayName', ['Axon Segment 5: ', interval1_range]);
        else
            plot(t_vec, Vm_trace, 'b-', 'HandleVisibility', 'off');
        end
    end

    % --- Interval 2: Axon Segment 10 ---
    x_positions2 = linspace(interval2(1), interval2(2), 4);
    spatial_indices2 = round(x_positions2 / dx) + 1;

    % Interval range for legend
    interval2_range = sprintf('%.2f-%.4f cm', interval2(1), interval2(2));

    for i = 1:length(spatial_indices2)
        xi = spatial_indices2(i);
        Vm_trace = data.Vm_all(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'r-', 'DisplayName', ['Axon Segment 10: ', interval2_range]);
        else
            plot(t_vec, Vm_trace, 'r-', 'HandleVisibility', 'off');
        end
    end

    % --- Interval 3: Axon Segment 15 ---
    x_positions3 = linspace(interval3(1), interval3(2), 4);
    spatial_indices3 = round(x_positions3 / dx) + 1;

    % Interval range for legend
    interval3_range = sprintf('%.3f-%.4f cm', interval3(1), interval3(2));

    for i = 1:length(spatial_indices3)
        xi = spatial_indices3(i);
        Vm_trace = data.Vm_all(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'k-', 'DisplayName', ['Axon Segment 15: ', interval3_range]);
        else
            plot(t_vec, Vm_trace, 'k-', 'HandleVisibility', 'off');
        end
    end

    hold off
    xlabel('Time (ms)')
    ylabel('$V_m$ in millivolts', 'Interpreter', 'latex')
    title('Axon segment of $V_m$ vs Time for Set (1) on DC', 'Interpreter', 'latex')
    legend('show', 'Location', 'northeast')

    T = data.T;
    xmin = 0;
    xmax = T;
    ymin = -70;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits
end

function plot_Vm_minus_Vmy_picture(data, time_shot)

    L = data.L;
    m = data.m;
    dt = data.dt;

    % plotting Voltage vs Axon length
    figure(1)
    t = linspace(0, L, m);
    plot(t, data.Vm_all(round(time_shot/dt),:), 'b-');
    hold on
    plot(t, data.Vm_minus_Vmy(round(time_shot/dt),:), 'r-');
    hold off
    
    xmin = 0;
    xmax = L;
    ymin = -80;
    ymax = 40;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size

    text(xmin+0.121, ymax-26, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_m$', '$V_m - V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel("Length of the axon in cm.", 'FontSize', 16)
    title('$V_m$ and $V_m - V_{my}$ for Set (2) on DC model', 'Interpreter', 'latex', 'FontSize', 13);

    
    % Plotting just Vmy space
    figure(2)
    t = linspace(0, L, m);
    plot(t, data.Vmy_all(round(time_shot/dt),:), 'Color', [0 0.5 0], 'LineStyle', '-');
    hold off

    xmin = 0;
    xmax = L;
    ymin = -20;
    ymax = 80;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size

    text(xmin+0.121, ymax-16, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$V_{my}$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel("Length of the axon in cm.", 'FontSize', 16)
    title('$V_{my}$ for Set (1) on SC model', 'Interpreter', 'latex', 'FontSize', 13);
    
    
    % Plotting V_m, V_m-V_my, and V_my all together
    figure(3)
    L = data.L;
    T = data.T;
    m = data.m;
    dt = data.dt;
    
    xmin = 0;
    xmax = T;
    ymin = -70;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits

    % plotting Voltage vs Axon length
    figure(3)
    t = linspace(0, L, m);
    plot(t, data.Vm_all(round(time_shot/dt),:), 'k-');
    hold on
    plot(t, data.Vm_minus_Vmy(round(time_shot/dt),:), 'r-');
    hold on
    plot(t, data.Vmy_all(round(time_shot/dt),:), 'b-');
    hold off

    xmin = 0;
    xmax = L;
    ymin = -70;
    ymax = 50;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size

    text(xmin+0.121, ymax-29, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_m$', '$V_m - V_{my}$', '$V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel("Length of the axon in cm.", 'FontSize', 16)
    title('$V_m$ and $V_m - V_{my}$ for Set (2) on DC model', 'Interpreter', 'latex', 'FontSize', 13);

end

% Animation that plots the voltage vs time
function plot_animation_voltage_vs_time(data, p)
    T = data.T;
    m = data.m;
    n = data.n;
    dx = data.dx;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n); 
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data, 'Vmy_all') && isfield(data, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    for j=1:length(voltages)
        figure(1);
        hold on;

        xmin = 0;
        xmax = T;
        ymin = -90;
        ymax = 90;

        if strcmp(voltages{j}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{j}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end

        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Time in milliseconds');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:m
            x = data.(voltages{j})(:,i); 

            % Plot the vector
            plot(t, x, 'b-');
            hold on

            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

% Animation that plots the voltage vs space
function plot_animation_voltage_vs_space(data, p)
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;
    
    % SPATIAL PROFILE %
    % x axis is the axon length
    t = linspace(0, L, m); 
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data, 'Vmy_all') && isfield(data, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    % plotting animation of Vm (or if SC or DC, Vm, Vmy and Vm-Vmy)
    for j=1:length(voltages)
        figure(1);
        hold on;

        xmin = 0;
        xmax = L;
        ymin = -90;
        ymax = 90;
        
        if strcmp(voltages{j}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{j}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end
            
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Axon length in cm');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:n
            x = data.(voltages{j})(i,:);
            
            % Plot the vector
            plot(t, x, 'b-');
            
            % legend('DC model: $R_{pa} = 1000 k\Omega cm$', 'Location', 'northeast', 'Interpreter', 'latex')
            text(xmin + 0.05, ymax + 0.8, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
    
end

% Animation that plots the probabilities vs time
function plot_animation_probabilities_vs_space(data, p)
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = 0;
    ymax = 1;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Axon length in cm.');
    ylabel('Probabilities of gates begin open.');
    
    % Loop through each vector and plot them one by one
    for i = 1:m
        x1 = data.N_all(i,:);  
        x2 = data.M_all(i,:);  
        x3 = data.H_all(i,:);  

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        p3 = plot(t, x3, 'g-');
        hold on
        
        % Add the legend (NOTE: legends are what slows down animation)
        legend([p1, p2, p3], 'K activation gate', 'Na activation gate', 'Na inactivation date', 'Location', 'northeast');

        text(xmin, ymax, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
        
        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

% Animation that plots the probabilities vs time
function plot_animation_probabilities_vs_time(data, p)
    T = data.T;
    m = data.m;
    n = data.n;
    dx = data.dx;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = 0;
    ymax = 1;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Time in milliseconds');
    ylabel('Probabilities of gates begin open.');
    
    % Loop through each vector and plot them one by one
    for i = 1:m
        x1 = data.N_all(:,i);  
        x2 = data.M_all(:,i);  
        x3 = data.H_all(:,i);  

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        p3 = plot(t, x3, 'g-');
        hold on
        
        % Add the legend (NOTE: legends are what slows down animation)
        legend([p1, p2, p3], 'K activation gate', 'Na activation gate', 'Na inactivation date', 'Location', 'northeast');

        text(xmin, ymax + 0.05, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');
        
        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end



% plotting multiple time and space shots of voltage for squid
function plot_time_and_space_shots(data, list_of_positions, list_of_times)

    L = data.L;
    T = data.T;
    m = data.m;
    n = data.n;
    dt = data.dt;
    dx = data.dx;
    
    % plotting Voltage vs Axon length
    figure(1)
    t1 = linspace(0, L, m);
    for i = 1:length(list_of_times)
        plot(t1, data.Vm_all(round(list_of_times{i}/dt),:))
        hold on
    end
    
    % describing plots using legends
    legendStrings1 = {};
    for i  = 1:length(list_of_times)
        legendStrings1{end+1} = sprintf('$V_m$ at t = %.2f ms', list_of_times{i});
    end
    legend(legendStrings1, 'Interpreter','latex', 'FontSize', 14)
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel("Length of the axon in cm.", 'FontSize', 16)
    title('SC Myelianted: $V_m$ vs Space', 'Interpreter', 'latex', 'FontSize', 13);
    
    xmin = 0;
    xmax = L;
    ymin = -65;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size

    % plotting Voltage vs Time
    figure(2)
    
    
    t2 = linspace(0, T, n); % FULL MATRIX
    for i = 1:length(list_of_positions)
        plot(t2, data.Vm_all(:,round(list_of_positions{i}/dx)))
        hold on
    end
    
    % describing plots using legends
    legendStrings2 = {};
    for i = 1:length(list_of_positions)
        legendStrings2{end+1} = sprintf('$V_m$ at x = %.4f cm', list_of_positions{i});
    end

    xmin = 0;
    xmax = T;
    ymin = -65;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size
    
    legend(legendStrings2, 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$V_m$ in millivolts', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel("Time (ms)", 'FontSize', 16)
    title('Temporal Profile of $V_m$ for Set (2) on DC model', 'Interpreter', 'latex', 'FontSize', 13);
    
    % plotting N, M, H probability vs time (at the first position list_of_positions(1))
    figure(3)
    plot(t2, data.N_all(:,round(list_of_positions{1}/dx)))
    hold on
    plot(t2, data.M_all(:,round(list_of_positions{1}/dx)))
    hold on
    plot(t2, data.H_all(:,round(list_of_positions{1}/dx)))
    legendStrings3 = {
        sprintf('n at x = %g cm', list_of_positions{1}), ...
        sprintf('m at x = %g cm', list_of_positions{1}), ...
        sprintf('h at x = %g cm', list_of_positions{1})};
    legend(legendStrings3, 'Interpreter','latex', 'FontSize', 14)
    ylabel("Probabilities of ion channels opening/closing.", 'FontSize', 13)
    xlabel("Time in milliseconds.", 'FontSize', 13)
end

% Plots animation of axon voltage in two different experiments
function plot_voltage_vs_time_comparison_animation(data_set, p)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    T = data_set{1}.T;
    m = data_set{1}.m;
    n = data_set{1}.n;
    dx = data_set{1}.dx;
    
    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n);
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data_set{1}, 'Vmy_all') && isfield(data_set{1}, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 90;
    
    for k=1:length(voltages)
        
        if strcmp(voltages{k}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{k}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end

        axis([xmin xmax ymin ymax]);  % Set axis limits
        set(gca, 'FontSize', 13); % Set axis size
        xlabel('Time in milliseconds', 'FontSize', 16);
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex', 'FontSize', 16)
        
        % Loop through each vector and plot them one by one
        for i = 1:m

            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(:,i));
                hold on
            end

            % legend('SiGe Tube params', 'Tube+Paralyne params', 'Location', 'northeast');
            legend('DC Set (1) r_{pn}=50', 'DC Set (1) r_{pn}=1', 'DC Set (2) r_{pn}=50', 'DC Set (2) r_{pn}=1', 'Location', 'northeast', 'FontSize', 14);
            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

function plot_voltage_vs_space_comparison_animation(data_set, p)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    L = data_set{1}.L;
    m = data_set{1}.m;
    n = data_set{1}.n;
    dt = data_set{1}.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 
       
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data_set{1}, 'Vmy_all') && isfield(data_set{1}, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    for k=1:length(voltages)
        
        if strcmp(voltages{k}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{k}, 'Vmy_all')
            voltage_name = 'V_{my}';
            ymin = -20
            ymax = 120
        else
            voltage_name = 'V_m - V_{my}';
        end
    
        axis([xmin xmax ymin ymax]);  % Set axis limits
        set(gca, 'FontSize', 13); % Set axis size
        xlabel('Length of axon in cm');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:n
            % for j = 1:length(data_set)
            %     plot(t, data_set{j}.(voltages{k})(i,:));
            %     hold on
            % end
            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(i,:));
                hold on
            end

            % Add the legend (NOTE: the legend is what is slowing down the animation)
            legend('DC', 'DC Scaled Condcutances', 'Location', 'northeast', 'FontSize', 14);
            % legend('DC: Myelinated', 'DC: Tube', 'DC: Tube+Paralyene', 'Location', 'northeast');
            % legend('SC: Tube', 'DC: Tube', 'Location', 'northeast');
            % legend('SC model: Cohen DC Params', 'DC model: $R_{pa}, R_{pn}$ given', 'DC model: $R_{pa}, R_{pn}$ computed', 'Location', 'northeast', 'Interpreter', 'latex');
            % legend('SC model: $R_i = 0.144 k\Omega cm$', 'DC model: $R_i = 0.712 k\Omega cm$', 'Location', 'northeast', 'Interpreter', 'latex')
            % legend('DC model: $R_{pa} = 0.01 k\Omega cm$', 'DC model: $R_{pa} = 0.10 k\Omega cm$', 'DC model: $R_{pa} = 0.20 k\Omega cm$', 'DC model: $R_{pa} = 0.30 k\Omega cm$', 'DC model: $R_{pa} = 0.40 k\Omega cm$', 'DC model: $R_{pa} = 0.50 k\Omega cm$', 'Location', 'northeast', 'Interpreter', 'latex')
            % legend('SC model w/ DC params', 'DC model: $R_{pa} = 0.5 k\Omega cm$', 'DC model: $R_{pa} = 5 k\Omega cm$', 'DC model: $R_{pa} = 200 k\Omega cm$', 'DC model: $R_{pa} = 1000 k\Omega cm$', 'Location', 'northeast', 'Interpreter', 'latex');
            % legend('SC: cell 6 params', 'SC: avg. params', 'SC: DC cell 6 params', 'DC: cell 6 params', 'DC: avg. params', 'Location', 'northeast', 'Interpreter', 'latex');
            text(xmin, ymax + 0.1, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

% Plots animation of axon voltage in two different experiments
function plot_voltage_vs_time_comparison(data_set, space_shot)
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    T = data_set{1}.T;
    dx = data_set{1}.dx;
    n = data_set{1}.n;
    
    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n);
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 90;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size
    text(xmin+19.8, ymax-40, sprintf('Position: %.5f cm', space_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    xlabel('Time in milliseconds', 'FontSize', 16);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)

    for j = 1:length(data_set)
        plot(t, data_set{j}.Vm_all(:,round(space_shot/dx)));
        hold on
    end

    % legend('SiGe Tube params', 'Tube+Paralyne params', 'Location', 'northeast');
    legend('DC Set (1) r_pn 50', 'DC Set (2) r_pn 50', 'Location', 'northeast', 'FontSize', 14);
  
end

function plot_voltage_vs_space_comparison(data_set, time_shot)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    L = data_set{1}.L;
    m = data_set{1}.m;
    dt = data_set{1}.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    axis([xmin xmax ymin ymax]);  % Set axis limits
    set(gca, 'FontSize', 13); % Set axis size
    text(xmin+19.8, ymax-40, sprintf('Position: %.3f cm', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    xlabel('Length of axon in cm', 'FontSize', 16);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 16)

        
    for j = 1:length(data_set)
        plot(t, data_set{j}.Vm_all(round(time_shot/dt),:));
        hold on
    end

    % Add the legend (NOTE: the legend is what is slowing down the animation)
    legend('DC', 'DC Scaled Condcutances', 'Location', 'northeast', 'FontSize', 14);
    
end





% plotting Vm vs Space and Vm-Vmy vs Space in the same plot (only for SC and DC)
function plot_Vm_and_Vm_minus_Vmy_vs_space(data, p)
    
    % defining spatial and time variables from one of the data sets
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % x axis is the axon time
    t = linspace(0, L, m); 
    
    figure(1);
    
    title('1) SC: cell 6 params Spatial Profile', 'Interpreter', 'latex')

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    % Loop through each vector and plot them one by one
    for i = 1:n
        cla;
        
        plot(t, data.Vm_all(i,:), 'b-');
        hold on
        plot(t, data.Vm_minus_Vmy(i,:), 'r-');
        hold off
        
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Length of axon in cm');
        ylabel('Voltage in millivolts.', 'Interpreter', 'latex')
        title('1) DC: cell 6 params Spatial Profile', 'Interpreter', 'latex');

        % Add the legend (NOTE: the legend is what is slowing down the animation)
        legend('$V_m$', '$V_m - V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex');
        text(xmin+0.01, ymax-10, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
    
        % Add a pause to create animation effect
        pause(p);

    end
end

function plot_voltage_vs_space_comparison_variable_dt(data_set, p)

    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    L = data_set{1}.L;
    m = data_set{1}.m;
    n = data_set{1}.n;
    dt = data_set{1}.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 
       
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data_set{1}, 'Vmy_all') && isfield(data_set{1}, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    for k=1:length(voltages)
        
        if strcmp(voltages{k}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{k}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end
    
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Length of axon in cm');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:n
            

            plot(t, data_set{1}.(voltages{k})(i,:), 'b-');
            hold on
            plot(t, data_set{2}.(voltages{k})(i*10,:), 'r--');
            hold on 

            % Add the legend (NOTE: the legend is what is slowing down the animation)
            legend('DC: Myelinated dt = 0.01', 'DC: Myelinated dt = 0.001', 'Location', 'northeast');
            text(xmin, ymax + 0.1, sprintf('Time: %.2f ms', round(i*dt, 2)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end
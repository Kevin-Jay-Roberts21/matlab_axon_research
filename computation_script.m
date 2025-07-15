% A script used to run computation scripts on cretain datas
% Kevin Roberts
% July 2025

clear all
close all
clc

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

DC_Cohen_set1_Cmy01 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy01.mat');
DC_Cohen_set1_Cmy02 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy02.mat');
DC_Cohen_set1_Cmy03 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy03.mat');
DC_Cohen_set1_Cmy04 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy04.mat');
DC_Cohen_set1_Cmy05 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy05.mat');
DC_Cohen_set1_Cmy06 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy06.mat');
DC_Cohen_set1_Cmy07 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy07.mat');
DC_Cohen_set1_Cmy08 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy08.mat');
DC_Cohen_set1_Cmy09 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy09.mat');
DC_Cohen_set1_Cmy10 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy10.mat');
DC_Cohen_set1_Cmy11 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy11.mat');
DC_Cohen_set1_Cmy12 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy12.mat');
DC_Cohen_set1_Cmy13 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy13.mat');
DC_Cohen_set1_Cmy14 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy14.mat');
DC_Cohen_set1_Cmy15 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy15.mat');
DC_Cohen_set1_Cmy16 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy16.mat');
DC_Cohen_set1_Cmy17 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy17.mat');
DC_Cohen_set1_Cmy18 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy18.mat');
DC_Cohen_set1_Cmy19 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy19.mat');
DC_Cohen_set1_Cmy20 = load('projects/axon_simulations/Paper_Parameter_Sensitivity/DC_Cohen_set1_Cmy20.mat');
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





% different spaces used for CV calculation
HH_space_vec = [4 5; 5 6; 6 7; 7 8;]; % all in (cm)
HH_space_vec1 = [2 2.5; 2.5 3; 3 3.5; 3.5 4];
SC_and_DC_space_vec = [0.06 0.08; 0.08 0.1; 0.1 0.12; 0.12 0.14;]; % all in (cm)
spaces = SC_and_DC_space_vec;

% data = {DC_Cohen_set1_rpn100};
% data = {DC_Cohen_set1_rpn50, DC_Cohen_set1_rpn200, DC_Cohen_set1_rpn400, DC_Cohen_set1_rpn600, DC_Cohen_set1_rpn800, DC_Cohen_set1_rpn1000, DC_Cohen_set1_rpn1200, DC_Cohen_set1_rpn1400, DC_Cohen_set1_rpn1600, DC_Cohen_set1_rpn1800, DC_Cohen_set1_rpn2000, DC_Cohen_set1_rpn2200, DC_Cohen_set1_rpn2400, DC_Cohen_set1_rpn2600, DC_Cohen_set1_rpn2800, DC_Cohen_set1_rpn3000};
% data = {DC_Cohen_set2_rpn50, DC_Cohen_set2_rpn200, DC_Cohen_set2_rpn400, DC_Cohen_set2_rpn600, DC_Cohen_set2_rpn800, DC_Cohen_set2_rpn1000, DC_Cohen_set2_rpn1200, DC_Cohen_set2_rpn1400, DC_Cohen_set2_rpn1600, DC_Cohen_set2_rpn1800, DC_Cohen_set2_rpn2000, DC_Cohen_set2_rpn2200, DC_Cohen_set2_rpn2400, DC_Cohen_set2_rpn2600, DC_Cohen_set2_rpn2800, DC_Cohen_set2_rpn3000};
data = {DC_Cohen_set1_Cmy01, DC_Cohen_set1_Cmy02, DC_Cohen_set1_Cmy03, DC_Cohen_set1_Cmy04, DC_Cohen_set1_Cmy05, DC_Cohen_set1_Cmy06, DC_Cohen_set1_Cmy07, DC_Cohen_set1_Cmy08, DC_Cohen_set1_Cmy09, DC_Cohen_set1_Cmy10, DC_Cohen_set1_Cmy11, DC_Cohen_set1_Cmy12, DC_Cohen_set1_Cmy13, DC_Cohen_set1_Cmy14, DC_Cohen_set1_Cmy15, DC_Cohen_set1_Cmy16, DC_Cohen_set1_Cmy17, DC_Cohen_set1_Cmy18, DC_Cohen_set1_Cmy19, DC_Cohen_set1_Cmy20};
% data = {DC_Cohen_set2_Cmy01, DC_Cohen_set2_Cmy02, DC_Cohen_set2_Cmy03, DC_Cohen_set2_Cmy04, DC_Cohen_set2_Cmy05, DC_Cohen_set2_Cmy06, DC_Cohen_set2_Cmy07, DC_Cohen_set2_Cmy08, DC_Cohen_set2_Cmy09, DC_Cohen_set2_Cmy10, DC_Cohen_set2_Cmy11, DC_Cohen_set2_Cmy12, DC_Cohen_set2_Cmy13, DC_Cohen_set2_Cmy14, DC_Cohen_set2_Cmy15, DC_Cohen_set2_Cmy16, DC_Cohen_set2_Cmy17, DC_Cohen_set2_Cmy18, DC_Cohen_set2_Cmy19, DC_Cohen_set2_Cmy20};
% data = {DC_Cohen_set1_Rmy50, DC_Cohen_set1_Rmy100, DC_Cohen_set1_Rmy200, DC_Cohen_set1_Rmy300, DC_Cohen_set1_Rmy400, DC_Cohen_set1_Rmy500, DC_Cohen_set1_Rmy600, DC_Cohen_set1_Rmy700, DC_Cohen_set1_Rmy800, DC_Cohen_set1_Rmy900, DC_Cohen_set1_Rmy1000};
% data = {DC_Cohen_set2_Rmy50, DC_Cohen_set2_Rmy100, DC_Cohen_set2_Rmy200, DC_Cohen_set2_Rmy300, DC_Cohen_set2_Rmy400, DC_Cohen_set2_Rmy500, DC_Cohen_set2_Rmy600, DC_Cohen_set2_Rmy700, DC_Cohen_set2_Rmy800, DC_Cohen_set2_Rmy900, DC_Cohen_set2_Rmy1000};


% printing the peak Vm of AP, and AP duration, and conduction velocity
AP_peaks = computing_AP_peak(data);
AP_durations = computing_AP_duration(data);
cvs = computing_conduction_velocity(spaces, data);

fprintf('Action Potential Voltage Peak: %.2f mV\n', computing_AP_peak(data));
fprintf('Action Potential Duration: %.2f ms\n', computing_AP_duration(data));
fprintf('Conduction Velocity: %.2f m/s\n', computing_conduction_velocity(spaces, data));
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
SC_Cohen_set1_new_conductances = load('projects/axon_simulations/Paper_SC_DC_Data/SC_Cohen_set1_new_conductances.mat');
% DC_Cohen_set1 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1.mat');
% DC_Cohen_set1_long = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_long.mat');
% DC_Cohen_set2 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set2.mat');
% DC_Cohen_set1_T33 = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_T33.mat');
% DC_Cohen_set1_T33_long = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_T33_long.mat');
DC_Cohen_set1_new_conductances = load('projects/axon_simulations/Paper_SC_DC_Data/DC_Cohen_set1_new_conductances.mat');


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



% different spaces used for CV calculation
HH_space_vec = [4 5; 5 6; 6 7; 7 8;]; % all in (cm)
HH_space_vec1 = [2 2.5; 2.5 3; 3 3.5; 3.5 4];
SC_and_DC_space_vec = [0.06 0.08; 0.08 0.1; 0.1 0.12; 0.12 0.14;]; % all in (cm)
spaces = SC_and_DC_space_vec;


data = {SC_Cohen_set1_new_conductances, DC_Cohen_set1_new_conductances};

% printing the peak Vm of AP, and AP duration, and conduction velocity
fprintf('Action Potential Voltage Peak: %.2f mV\n', computing_AP_peak(data));
fprintf('Action Potential Duration: %.2f ms\n', computing_AP_duration(data));
fprintf('Conduction Velocity: %.2f m/s\n', computing_conduction_velocity(spaces, data));
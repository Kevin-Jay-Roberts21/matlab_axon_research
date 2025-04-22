% A MATLAB script to compute the conduction velocity of the voltage
% simulations. More accurate results for longer simulations.
% Kevin Roberts
% February 2025

clear all
close all
clc

% Collecting temp data for HH
% HH_temp_base = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_6.3.mat");
% HH_temp_8 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_8.mat");
% HH_temp_10 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_10.mat");
% HH_temp_12 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_12.mat");
% HH_temp_14 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_14.mat");
% HH_temp_16 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_16.mat");
% HH_temp_18 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_18.mat");
% HH_temp_20 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_20.mat");
% HH_temp_22 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_22.mat");
% HH_temp_24 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_24.mat");
% HH_temp_26 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_26.mat");
% HH_temp_28 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_28.mat");
% HH_temp_30 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_30.mat");
% HH_temp_31 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_31.mat");
% HH_temp_32 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_32.mat");
% HH_temp_33 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_33.mat");
% HH_temp_34 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_34.mat");
% HH_temp_35 = load("projects/axon_simulations/HH_temp_data_stim_20/HH_Temp_35.mat");

% Addition temp data for HH
% HH_data_Temp_base = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_6.3.mat');
% HH_data_Temp_7 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_7.mat');
% HH_data_Temp_8 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_8.mat');
% HH_data_Temp_9 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_9.mat');
% HH_data_Temp_10 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_10.mat');
% HH_data_Temp_11 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_11.mat');
% HH_data_Temp_12 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_12.mat');
% HH_data_Temp_13 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_13.mat');
% HH_data_Temp_14 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_14.mat');
% HH_data_Temp_15 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_15.mat');
% HH_data_Temp_16 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_16.mat');
% HH_data_Temp_17 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_17.mat');
% HH_data_Temp_18 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_18.mat');
% HH_data_Temp_19 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_19.mat');
% HH_data_Temp_20 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_20.mat');
% HH_data_Temp_21 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_21.mat');
% HH_data_Temp_22 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_22.mat');
% HH_data_Temp_23 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_23.mat');
% HH_data_Temp_24 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_24.mat');
% HH_data_Temp_25 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_25.mat');
% HH_data_Temp_26 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_26.mat');
% HH_data_Temp_27 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_27.mat');
% HH_data_Temp_28 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_28.mat');
% HH_data_Temp_29 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_29.mat');
% HH_data_Temp_30 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_30.mat');
% HH_data_Temp_31 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_31.mat');
% HH_data_Temp_32 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_32.mat');
% HH_data_Temp_33 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_33.mat');
% HH_data_Temp_34 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_34.mat');
% HH_data_Temp_35 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_35.mat');
% 
% SC temp simulations on DC cell6 params 
% SC_temp_20 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_base.mat');
% SC_temp_22 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_22.mat');
% SC_temp_24 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_24.mat');
% SC_temp_26 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_26.mat');
% SC_temp_28 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_28.mat');
% SC_temp_30 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_30.mat');
% SC_temp_32 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_32.mat');
% SC_temp_34 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_34.mat');
% SC_temp_36 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_36.mat');
% SC_temp_38 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_38.mat');
% SC_temp_40 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_40.mat');
% SC_temp_42 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_42.mat');
% SC_temp_44 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_44.mat');
% SC_temp_46 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_46.mat');
% SC_temp_48 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_48.mat');
% SC_temp_50 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_50.mat');
% SC_temp_51 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_51.mat');
% SC_temp_52 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_52.mat');
% SC_temp_53 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_53.mat');
% SC_temp_54 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_54.mat');
% SC_temp_55 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_55.mat');
% SC_temp_56 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_56.mat');
% SC_temp_57 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_57.mat');
% SC_temp_58 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_58.mat');

% SC temp simulations on DC cell6 params 
% SC_temp_20 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_base.mat');
% SC_temp_22 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_22.mat');
% SC_temp_24 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_24.mat');
% SC_temp_26 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_26.mat');
% SC_temp_28 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_28.mat');
% SC_temp_30 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_30.mat');
% SC_temp_32 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_32.mat');
% SC_temp_34 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_34.mat');
% SC_temp_36 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_36.mat');
% SC_temp_38 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_38.mat');
% SC_temp_40 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_40.mat');
% SC_temp_42 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_42.mat');
% SC_temp_44 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_44.mat');
% SC_temp_46 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_46.mat');
% SC_temp_48 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_48.mat');
% SC_temp_50 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_50.mat');
% SC_temp_51 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_51.mat');
% SC_temp_52 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_52.mat');
% SC_temp_53 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_53.mat');
% SC_temp_54 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_54.mat');
% SC_temp_55 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_55.mat');
% SC_temp_56 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_56.mat');
% SC_temp_57 = load('projects/axon_simulations/SC_temp_data_on_DCcell6_params/SC_temp_57.mat');

% DC temp simulations on DC cell6 params 
% DC_temp_20 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_base.mat');
% DC_temp_21 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_21.mat');
% DC_temp_22 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_22.mat');
% DC_temp_23 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_23.mat');
% DC_temp_24 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_24.mat');
% DC_temp_25 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_25.mat');
% DC_temp_26 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_26.mat');
% DC_temp_27 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_27.mat');
% DC_temp_28 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_28.mat');
% DC_temp_29 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_29.mat');
% DC_temp_30 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_30.mat');
% DC_temp_31 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_31.mat');
% DC_temp_32 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_32.mat');
% DC_temp_33 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_33.mat');

% SC temp simulations on DC avg params 
% SC_temp_20 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_base.mat');
% SC_temp_22 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_22.mat');
% SC_temp_24 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_24.mat');
% SC_temp_26 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_26.mat');
% SC_temp_28 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_28.mat');
% SC_temp_30 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_30.mat');
% SC_temp_32 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_32.mat');
% SC_temp_34 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_34.mat');
% SC_temp_36 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_36.mat');
% SC_temp_38 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_38.mat');
% SC_temp_40 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_40.mat');
% SC_temp_42 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_42.mat');
% SC_temp_44 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_44.mat');
% SC_temp_46 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_46.mat');
% SC_temp_48 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_48.mat');
% SC_temp_50 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_50.mat');
% SC_temp_51 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_51.mat');
% SC_temp_52 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_52.mat');
% SC_temp_53 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_53.mat');
% SC_temp_54 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_54.mat');
% SC_temp_55 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_55.mat');
% SC_temp_56 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_56.mat');
% SC_temp_57 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_57.mat');

% DC temp simulations on DC avg params 
% DC_temp_20 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_base.mat');
% DC_temp_21 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_21.mat');
% DC_temp_22 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_22.mat');
% DC_temp_23 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_23.mat');

% Testing different stimulus values for SC
% SC_Cohen_DC_cell6_stim206 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim206.mat');
% SC_Cohen_DC_cell6_stim306 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim306.mat');
% SC_Cohen_DC_cell6_stim406 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim406.mat');
% SC_Cohen_DC_cell6_stim506 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim506.mat');
% SC_Cohen_DC_cell6_stim606 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim606.mat');
% SC_Cohen_DC_cell6_stim706 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim706.mat');
% SC_Cohen_DC_cell6_stim806 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim806.mat');
% SC_Cohen_DC_cell6_stim906 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim906.mat');
% SC_Cohen_DC_cell6_stim1006 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1006.mat');
% SC_Cohen_DC_cell6_stim1106 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1106.mat');
% SC_Cohen_DC_cell6_stim1206 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1206.mat');
% SC_Cohen_DC_cell6_stim1306 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1306.mat');
% SC_Cohen_DC_cell6_stim1406 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1406.mat');
% SC_Cohen_DC_cell6_stim1506 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1506.mat');
% SC_Cohen_DC_cell6_stim1606 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1606.mat');
% SC_Cohen_DC_cell6_stim1706 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1706.mat');
% SC_Cohen_DC_cell6_stim1806 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1806.mat');
% SC_Cohen_DC_cell6_stim1906 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim1906.mat');
% SC_Cohen_DC_cell6_stim2006 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim2006.mat');
% SC_Cohen_DC_cell6_stim3006 = load('projects/axon_simulations/Stimulus Increase/SC_Cohen_DC_cell6_stim3006.mat');

% Testing different stimulus values for DC
% DC_Cohen_DC_cell6_stim846 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim846.mat');
% DC_Cohen_DC_cell6_stim946 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim946.mat');
% DC_Cohen_DC_cell6_stim1046 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1046.mat');
% DC_Cohen_DC_cell6_stim1146 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1146.mat');
% DC_Cohen_DC_cell6_stim1246 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1246.mat');
% DC_Cohen_DC_cell6_stim1346 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1346.mat');
% DC_Cohen_DC_cell6_stim1446 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1446.mat');
% DC_Cohen_DC_cell6_stim1546 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1546.mat');
% DC_Cohen_DC_cell6_stim1646 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1646.mat');
% DC_Cohen_DC_cell6_stim1746 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1746.mat');
% DC_Cohen_DC_cell6_stim1846 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1846.mat');
% DC_Cohen_DC_cell6_stim1946 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim1946.mat');
% DC_Cohen_DC_cell6_stim2046 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim2046.mat');
% DC_Cohen_DC_cell6_stim3046 = load('projects/axon_simulations/Stimulus Increase/DC_Cohen_DC_cell6_stim3046.mat');

% Huang Parameters
SC_Huang_Myelinated_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_Myelinated_set1.mat');
SC_Huang_Tube_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_Tube_set1.mat');
SC_Huang_TubeParalyne_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_TubeParalyne_set1.mat');
% DC_Huang_Myelinated_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_set1.mat');
% DC_Huang_Tube_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Tube_set1.mat');
% DC_Huang_TubeParalyne_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_TubeParalyne_set1.mat');


% Cohen parameters
% SC_Cohen_DC_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_params.mat');
% SC_Cohen_DC_avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_avg_params.mat');
% DC_Cohen_DC_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params.mat');
% DC_Cohen_DC_avg_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_avg_params.mat');

% NOTE: Time and space choices will vary depending on HH, SC, and DC models

% creating a list of spaces (x1 and x2) to compute the cv, we will then
% take an average of these different spaces to get a more accurate cv
% NOTE: these may differ for different length axons
HH_space_vec = [4 5; 5 6; 6 7; 7 8;]; % all in (cm)
HH_space_vec1 = [2 2.5; 2.5 3; 3 3.5; 3.5 4];
SC_and_DC_space_vec = [0.06 0.08; 0.08 0.1; 0.1 0.12; 0.12 0.14;]; % all in (cm)

% defining multiple datasets
% set_of_data1 = {HH_temp_base, HH_temp_8, HH_temp_10, HH_temp_12, HH_temp_14, HH_temp_16, HH_temp_18, HH_temp_20, HH_temp_22, HH_temp_24, HH_temp_26, HH_temp_28, HH_temp_30, HH_temp_31, HH_temp_32};
% set_of_data2 = {SC_temp_20, SC_temp_22, SC_temp_24, SC_temp_26, SC_temp_28, SC_temp_30, SC_temp_32, SC_temp_34, SC_temp_36, SC_temp_38, SC_temp_40, SC_temp_42, SC_temp_44, SC_temp_46, SC_temp_48, SC_temp_50, SC_temp_52, SC_temp_54, SC_temp_55};
% set_of_data3 = {DC_temp_20, DC_temp_21, DC_temp_22, DC_temp_23};
set_of_data4 = {SC_Huang_Myelinated_set1, SC_Huang_Tube_set1, SC_Huang_TubeParalyne_set1};
% set_of_data5 = {SC_Cohen_DC_cell6_params, SC_Cohen_DC_avg_params, DC_Cohen_DC_cell6_params, DC_Cohen_DC_avg_params};
% set_of_data6 = {SC_temp_20, SC_temp_22, SC_temp_24, SC_temp_26, SC_temp_28, SC_temp_30, SC_temp_32, SC_temp_34, SC_temp_36, SC_temp_38, SC_temp_40, SC_temp_42, SC_temp_44, SC_temp_46, SC_temp_48, SC_temp_50, SC_temp_51, SC_temp_52, SC_temp_53, SC_temp_54, SC_temp_55, SC_temp_56, SC_temp_57};
% set_of_data7 = {DC_temp_20, DC_temp_21, DC_temp_22, DC_temp_23, DC_temp_24, DC_temp_25, DC_temp_26, DC_temp_27, DC_temp_28, DC_temp_29, DC_temp_30, DC_temp_31, DC_temp_32, DC_temp_33};
% picking the dataset to compute cv
data = set_of_data4;

% computing conduction velocity
cv = calculate_cv(SC_and_DC_space_vec, data) % computed cv's (in m/s)


% conduction velocity function
function final_cvs = calculate_cv(spaces, data)
    
    sum_of_cvs = zeros(1, length(data));

    for j = 1:size(spaces, 1)
        
        list_of_cv = zeros(1, length(data)); % Preallocate for efficiency
        
        for i = 1:length(data)
            
            x1 = spaces(j, 1);
            x2 = spaces(j, 2);
            
            % need to find the index of x1 and x2
            index_x1 = x1/data{i}.dx;
            index_x2 = x2/data{i}.dx;
            
            % identifying the space vectors at index_t1 and index_t2
            vec_t1 = data{i}.Vm_all(:,index_x1);
            vec_t2 = data{i}.Vm_all(:,index_x2);
        
            % computing the positions of where the max voltage is at t1 and t2:
            [Vm_at_x1, index_t1] = max(vec_t1); % Time index where voltage peaks at x1
            [Vm_at_x2, index_t2] = max(vec_t2); % Time index where voltage peaks at x2
            
            % Note that x1 and x2 are just the indices. We need to find their actual
            % spatial position in cm. This is done by identifiying the mesh from the
            % data
            t1 = index_t1*data{i}.dt; % (in cm)
            t2 = index_t2*data{i}.dt; % (in cm)
            
            % finally, compute the conduction velocity
            cv = (x2 - x1)/(t2 - t1) * 10; % *10 to convert to m/s
    
            list_of_cv(i) = cv;
        end
        
        sum_of_cvs = sum_of_cvs + list_of_cv;
    end

    % now dividing sum_of_cvs by the number of cv we computed
    final_cvs = sum_of_cvs/size(spaces, 1);

end




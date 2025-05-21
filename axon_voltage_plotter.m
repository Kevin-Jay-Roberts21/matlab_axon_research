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

% HH_1 = load('projects/axon_simulations/HH_basic/HH_1.mat');
% HH_2 = load('projects/axon_simulations/HH_basic/HH_2.mat');
% HH_3 = load('projects/axon_simulations/HH_basic/HH_3.mat');


% Huang parameters
% SC_Huang_Myelinated_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_Myelinated_set1.mat');
% SC_Huang_Tube_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_Tube_set1.mat');
% SC_Huang_TubeParalyene_set1 = load('projects/axon_simulations/Huang_simulations/SC_Huang_TubeParalyene_set1.mat');
% SC_Huang_Myelinated_set1_long = load('projects/axon_simulations/Huang_simulations/SC_Huang_Myelinated_set1_long.mat');
% SC_Huang_Tube_set1_long = load('projects/axon_simulations/Huang_simulations/SC_Huang_Tube_set1_long.mat');
% SC_Huang_TubeParalyene_set1_long = load('projects/axon_simulations/Huang_simulations/SC_Huang_TubeParalyene_set1_long.mat');
% DC_Huang_Myelinated_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_set1.mat');
% DC_Huang_Tube_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Tube_set1.mat');
% DC_Huang_TubeParalyne_set1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_TubeParalyne_set1.mat');
% DC_Huang_Myelianted_dt_01 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_dt_0.01.mat');
% DC_Huang_Myelianted_dt_001 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_dt_0.001.mat');
DC_Huang_Myelinated = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated.mat');
DC_Huang_Myelinated_test1 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_test1.mat');
DC_Huang_Myelinated_test2 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_test2.mat');
DC_Huang_Myelinated_test3 = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated_test3.mat');
SC_Huang_Myelianted_test3 = load('projects/axon_simulations/Huang_simulations/SC_Huang_Myelinated_test3.mat');

% Increasing R_pa
% SC_Cohen_DC_params = load('SC_model_with_DC_Cohen_params.mat');
% DC_Cohen_DC_params_Rpa_05 = load('DC_model_with_DC_Cohen_params_Rpa_0.5.mat');
% DC_Cohen_DC_params_Rpa_5 = load('DC_model_with_DC_Cohen_params_Rpa_5.mat');
% DC_Cohen_DC_params_Rpa_200 = load('DC_model_with_DC_Cohen_params_Rpa_200.mat');
% DC_Cohen_DC_params_Rpa_1000 = load('DC_model_with_DC_Cohen_params_Rpa_1000.mat');
% DC_Cohen_DC_params_Rpa_2000 = load('DC_model_with_DC_Cohen_params_Rpa_2000.mat');

% Cohen data
% SC_Cohen_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_SC_cell6_params.mat');
% SC_Cohen_avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_SC_avg_params.mat');
% SC_Cohen_DC_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_params.mat');
% SC_Cohen_DC_cell6_params_long = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_params_long.mat');
% SC_Cohen_DC_avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_avg_params.mat');
% SC_Cohen_DC_cell6_temp_25 = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_temp_25.mat');
% SC_Cohen_DC_cell6_temp_33_long = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_temp33_long.mat');
% SC_Cohen_DC_cell6_temp_56_long = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_cell6_temp56_long.mat');
% DC_Cohen_DC_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params.mat');
% DC_Cohen_DC_cell6_params_long = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params_long.mat');
% DC_Cohen_DC_cell5_params_shifted_stimulus = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell5_params.mat');
% DC_Cohen_DC_cell6_params_shifted_stimulus = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params_shifted_stimulus.mat');
% DC_Cohen_DC_cell6_params_shifted_stimulus_2 = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params_shifted_stimulus_2.mat');
% DC_Cohen_DC_cell6_params_shifted_stimulus_3 = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_params_shifted_stimulus_3.mat');
% DC_Cohen_DC_avg_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_avg_params.mat');
% DC_Cohen_DC_cell6_temp_25 = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_temp_25.mat');
% DC_Cohen_DC_cell6_temp_33_long = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_DC_cell6_temp33_long.mat');
% DC_Cohen_avg_r_pa1000fold = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_avg_r_pa1000fold.mat');
% DC_Cohen_avg_stim_increase = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_avg_stim_increase.mat');

% increaing R_pa
% DC_data_Rpa_01 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.01.mat');
% DC_data_Rpa_02 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.02.mat');
% DC_data_Rpa_03 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.03.mat');
% DC_data_Rpa_04 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.04.mat');
% DC_data_Rpa_05 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.05.mat');
% DC_data_Rpa_06 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.06.mat');
% DC_data_Rpa_07 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.07.mat');
% DC_data_Rpa_08 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.08.mat');
% DC_data_Rpa_09 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.09.mat');
% DC_data_Rpa_10 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.1.mat');
% DC_data_Rpa_15 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.15.mat');
% DC_data_Rpa_20 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.2.mat');
% DC_data_Rpa_25 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.25.mat');
% DC_data_Rpa_30 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.3.mat');
% DC_data_Rpa_35 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.35.mat');
% DC_data_Rpa_40 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.4.mat');
% DC_data_Rpa_45 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.45.mat');
% DC_data_Rpa_50 = load('projects/axon_simulations/DC_R_pa_data/DC_model_Rpa_0.5.mat');

% Temp Analysis data 1
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

% Temp Analysis data 2
% HH_data_Temp_base = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_6.3.mat');
% HH_data_Temp_7 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_7.mat');
% HH_data_Temp_8 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_8.mat');
% HH_data_Temp_9 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_9.mat');
% HH_data_Temp_10 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_10.mat');
% HH_data_Temp_11 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_11.mat');
% HH_data_Temp_12 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_12.mat');
% HH_data_Temp_13 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_13.mat');
% HH_data_Temp_14 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_14.mat');
% HH_data_Temp_15 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_15.mat');
% HH_data_Temp_16 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_16.mat');
% HH_data_Temp_17 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_17.mat');
% HH_data_Temp_18 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_18.mat');
% HH_data_Temp_19 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_19.mat');
% HH_data_Temp_20 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_20.mat');
% HH_data_Temp_21 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_21.mat');
% HH_data_Temp_22 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_22.mat');
% HH_data_Temp_23 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_23.mat');
% HH_data_Temp_24 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_24.mat');
% HH_data_Temp_25 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_25.mat');
% HH_data_Temp_26 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_26.mat');
% HH_data_Temp_27 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_27.mat');
% HH_data_Temp_28 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_28.mat');
% HH_data_Temp_29 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_29.mat');
% HH_data_Temp_30 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_30.mat');
% HH_data_Temp_31 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_31.mat');
% HH_data_Temp_32 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_32.mat');
% HH_data_Temp_33 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_33.mat');
% HH_data_Temp_34 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_34.mat');
% HH_data_Temp_35 = load('projects/axon_simulations/HH_temp_data_stim_12/HH_data_Temp_35.mat');

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
% DC_temp_34 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_34.mat');
% DC_temp_35 = load('projects/axon_simulations/DC_temp_data_on_DCcell6_params/DC_temp_35.mat');

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
% SC_temp_58 = load('projects/axon_simulations/SC_temp_data_on_DCavg_params/SC_temp_58.mat');

% DC temp simulations on DC avg params 
% DC_temp_20 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_base.mat');
% DC_temp_21 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_21.mat');
% DC_temp_22 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_22.mat');
% DC_temp_23 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_23.mat');
% DC_temp_24 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_24.mat');
% DC_temp_25 = load('projects/axon_simulations/DC_temp_data_on_DCavg_params/DC_temp_25.mat');


% picking interval
interval1 = [0.0400 0.0485]; % interval is in cm
interval2 = [0.0800 0.0885]; % interval is in cm
interval3 = [0.1200 0.1285];
time_shot = 1.1; % in ms

% picking time shots
% time1 = 5.1; % in ms
% time2 = 9.5; % in ms
% time3 = 10; % in ms
time1 = 2; % in ms
time2 = 3.3; % in ms
time3 = 4; % in ms
list_of_times = {time1, time2, time3};

% picking space shots
% position1 = 1; % in cm
% position2 = 3; % in cm
% position3 = 5; % in cm

position1 = 0.0800; % in cm
position2 = 0.0828; % in cm
position3 = 0.0857; % in cm
position4 = 0.0885; % in cm

% position1 = 0.0005; % in cm
% position2 = 0.06; % in cm
% position4 = 0.16; % in cm

list_of_positions = {position1, position2, position4};


% picking pause (this controls the speed of the animation, the pause variable 
% is in seconds. It is how many seconds each time (or space) shot will be
% paused at). NOTE: the legends are what slows down the animations, may
% condsider getting rid of or adding them in certain cases.
p = 0.001;

% creating a set of data from multiple experiments used to plot animation
% (the first element in the set_of_data is darkred, then the proceeding elements get
% brighter and brighter until the last element which is the brightest red)
% set_of_data1 = {HH_data_Temp_base, HH_data_Temp_30, HH_data_Temp_35};
% set_of_data2 = {HH_data_Temp_26, HH_data_Temp_27, HH_data_Temp_28, HH_data_Temp_29, HH_data_Temp_30, HH_data_Temp_31, HH_data_Temp_32};
% set_of_data3 = {SC_data_from_saltcond2023_code, DC_data_from_saltcond2023_code};
% set_of_data4 = {SC_data_Temp_20, SC_data_Temp_40, SC_data_Temp_52};
% set_of_data5 = {DC_data_Temp_20, SC_data_Temp_40, DC_data_Temp_52};
% set_of_data6 = {DC_Huang_Myelinated, DC_Huang_Tube_params, DC_Huang_TubeParalyne_params};
% set_of_data7 = {SC_Cohen_DC_params, DC_Cohen_DC_params_Rpa_05, DC_Cohen_DC_params_Rpa_5, DC_Cohen_DC_params_Rpa_200, DC_Cohen_DC_params_Rpa_1000};
% set_of_data8 = {SC_Ri_144, DC_model_SC_params_w_Ri_144};
% set_of_data9 = {DC_data_Rpa_01, DC_data_Rpa_05, DC_data_Rpa_10, DC_data_Rpa_20, DC_data_Rpa_30, DC_data_Rpa_40, DC_data_Rpa_50};
% set_of_data10 = {SC_Cohen_cell6_params, SC_Cohen_Avg_params, SC_Cohen_DC_cell6_params, DC_Cohen_cell6_params, DC_Cohen_Avg_params};
% set_of_data11 = {HH_temp_base, HH_temp_8, HH_temp_10, HH_temp_12, HH_temp_14, HH_temp_16, HH_temp_18, HH_temp_20, HH_temp_22, HH_temp_24, HH_temp_26, HH_temp_28, HH_temp_30, HH_temp_31, HH_temp_32, HH_temp_33, HH_temp_34, HH_temp_35};
% set_of_data12 = {SC_temp_52, SC_temp_53, SC_temp_54, SC_temp_55, SC_temp_56, SC_temp_57, SC_temp_58};
% set_of_data13 = {DC_temp_30, DC_temp_31, DC_temp_32, DC_temp_33, DC_temp_34, DC_temp_35};
set_of_data14 = {DC_Huang_Myelinated_test1, DC_Huang_Myelinated_test2, DC_Huang_Myelinated_test3, DC_Huang_Myelinated};

% data = DC_Cohen_DC_cell5_params_shifted_stimulus;
% data = SC_temp_58;

% plot_zoomed_in_region_w_AP_at_spaces(data, time_shot, interval1, interval2, interval3);
% plot_Vm_minus_Vmy_picture(data, time_shot);
% plot_animation_voltage_vs_time(data, p);
% plot_animation_voltage_vs_space(data, p);
% plot_animation_probabilities_vs_time(HH_data_Temp_33, p);
% plot_animation_probabilities_vs_space(HH_data_Temp_32, p);
% plot_time_and_space_shots(data, list_of_positions, list_of_times);
plot_voltage_vs_time_comparison(set_of_data14, p);
plot_voltage_vs_space_comparison(set_of_data14, p);
% plot_Vm_and_Vm_minus_Vmy_vs_space(data, p);
% plot_voltage_vs_space_comparison_variable_dt(set_of_data14, p);


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
    ymin = -70;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits

    text(xmin+0.10, ymax-30, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_m$', '$V_m - V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 13)
    xlabel("Length of the axon in cm.", 'FontSize', 13)
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

    text(xmin+0.12, ymax-30, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$V_{my}$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 13)
    xlabel("Length of the axon in cm.", 'FontSize', 13)
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

    text(xmin+0.11, ymax-30, sprintf('Time: %.3f ms', time_shot), 'FontSize', 11, 'BackgroundColor', 'w');
    legend('$V_m$', '$V_m - V_{my}$', '$V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 13)
    xlabel("Length of the axon in cm.", 'FontSize', 13)
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
    legend(legendStrings1, 'Interpreter','latex', 'FontSize', 12)
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex', 'FontSize', 13)
    xlabel("Length of the axon in cm.", 'FontSize', 13)
    title('SC Myelianted: $V_m$ vs Space', 'Interpreter', 'latex', 'FontSize', 13);
    
    xmin = 0;
    xmax = L;
    ymin = -70;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits

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
    ymin = -70;
    ymax = 30;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    legend(legendStrings2, 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$V_m$ (mV)', 'Interpreter', 'latex', 'FontSize', 13)
    xlabel("Time (ms)", 'FontSize', 13)
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
    legend(legendStrings3, 'Interpreter','latex', 'FontSize', 12)
    ylabel("Probabilities of ion channels opening/closing.", 'FontSize', 13)
    xlabel("Time in milliseconds.", 'FontSize', 13)
end

% Plots animation of axon voltage in two different experiments
function plot_voltage_vs_time_comparison(data_set, p)
    
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
        xlabel('Time in milliseconds');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')
        
        % Loop through each vector and plot them one by one
        for i = 1:m

            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(:,i));
                hold on
            end

            % legend('SiGe Tube params', 'Tube+Paralyne params', 'Location', 'northeast');
            legend('DC: Myelinated test1', 'DC: Myelinated test2', 'DC: Myelinated test3', 'DC: Myelinated', 'Location', 'northeast');
            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

function plot_voltage_vs_space_comparison(data_set, p)
    
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
            % for j = 1:length(data_set)
            %     plot(t, data_set{j}.(voltages{k})(i,:));
            %     hold on
            % end
            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(i,:));
                hold on
            end

            % Add the legend (NOTE: the legend is what is slowing down the animation)
            legend('DC: Myelinated test1', 'DC: Myelinated test2', 'DC: Myelinated test3', 'DC: Myelinated', 'Location', 'northeast');
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
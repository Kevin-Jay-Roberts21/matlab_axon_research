% Axon voltage plotter. This code is meant to plot data from the following
% files: 
% squid_giant_axon.m 
% myelianted_HH_axon.m,
% dimensionless_squid_giant_axon.m 
% sc_code_implicit_v1.m
% sc_code_implicit_v2.m
% sc_code_implicit_v3.m
% dc_code_implicit.m
% 
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

% getting the saved data
% HH_data = load('HH_data.mat');
HH_1 = load('HH_1.mat');
HH_2 = load('HH_2.mat');
HH_3 = load('HH_3.mat');
% HH_data2 = load('HH_data2.mat');
% SC_data = load('SC_data.mat');
% DC_data = load('DC_data.mat');

% SC_data1 = load('projects/SC_data1.mat');
% SC_data2 = load('projects/SC_data2.mat');
% SC_data3 = load('projects/SC_data3.mat');
% SC_data_Sv_200 = load('SC_data_Sv_200.mat');
% DC_data_Sv_200 = load('DC_data_Sv_200.mat');
% SC_data_from_saltcond2023_code = load('SC_data_from_saltcond2023_code.mat');
% DC_data_from_saltcond2023_code = load('DC_data_from_saltcond2023_code.mat');

% DC_Huang_Myelinated = load('projects/axon_simulations/Huang_simulations/DC_Huang_Myelinated.mat');
% DC_Huang_Tube_params = load('projects/axon_simulations/Huang_simulations/DC_Huang_Tube_params.mat');
% DC_Huang_TubeParalyne_params = load('projects/axon_simulations/Huang_simulations/DC_Huang_TubeParalyne_params.mat');

% SC_Ri_144 = load('SC_Ri_0.144.mat');

% SC_Cohen_DC_params = load('SC_model_with_DC_Cohen_params.mat');
% DC_Cohen_DC_params_Rpa_05 = load('DC_model_with_DC_Cohen_params_Rpa_0.5.mat');
% DC_Cohen_DC_params_Rpa_5 = load('DC_model_with_DC_Cohen_params_Rpa_5.mat');
% DC_Cohen_DC_params_Rpa_200 = load('DC_model_with_DC_Cohen_params_Rpa_200.mat');
% DC_Cohen_DC_params_Rpa_1000 = load('DC_model_with_DC_Cohen_params_Rpa_1000.mat');
% DC_Cohen_DC_params_Rpa_2000 = load('DC_model_with_DC_Cohen_params_Rpa_2000.mat');

% SC_Cohen_Optimized_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_Optimized_params.mat');
% SC_Cohen_Avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_Avg_params.mat');
% SC_Cohen_DC_Optimized_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_Optimized_params.mat');
% DC_Cohen_Optimized_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_Optimized_params.mat');
% DC_Cohen_Avg_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_Avg_params.mat');

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
% 
% SC_temp_35 = load('projects/axon_simulations/SC_temp_data/SC_temp_base.mat');
% SC_temp_37 = load('projects/axon_simulations/SC_temp_data/SC_temp_37.mat');
% SC_temp_39 = load('projects/axon_simulations/SC_temp_data/SC_temp_39.mat');
% SC_temp_41 = load('projects/axon_simulations/SC_temp_data/SC_temp_41.mat');
% SC_temp_43 = load('projects/axon_simulations/SC_temp_data/SC_temp_43.mat');
% SC_temp_45 = load('projects/axon_simulations/SC_temp_data/SC_temp_45.mat');
% SC_temp_47 = load('projects/axon_simulations/SC_temp_data/SC_temp_47.mat');
% SC_temp_49 = load('projects/axon_simulations/SC_temp_data/SC_temp_49.mat');
% SC_temp_51 = load('projects/axon_simulations/SC_temp_data/SC_temp_51.mat');
% SC_temp_53 = load('projects/axon_simulations/SC_temp_data/SC_temp_53.mat');
% SC_temp_55 = load('projects/axon_simulations/SC_temp_data/SC_temp_55.mat');
% SC_temp_57 = load('projects/axon_simulations/SC_temp_data/SC_temp_57.mat');
% SC_temp_59 = load('projects/axon_simulations/SC_temp_data/SC_temp_59.mat');
% SC_temp_61 = load('projects/axon_simulations/SC_temp_data/SC_temp_61.mat');
% SC_temp_62 = load('projects/axon_simulations/SC_temp_data/SC_temp_62.mat');
% SC_temp_63 = load('projects/axon_simulations/SC_temp_data/SC_temp_63.mat');
% SC_temp_64 = load('projects/axon_simulations/SC_temp_data/SC_temp_64.mat');
% SC_temp_65 = load('projects/axon_simulations/SC_temp_data/SC_temp_65.mat');
% SC_temp_66 = load('projects/axon_simulations/SC_temp_data/SC_temp_66.mat');
% SC_temp_67 = load('projects/axon_simulations/SC_temp_data/SC_temp_67.mat');
% SC_temp_68 = load('projects/axon_simulations/SC_temp_data/SC_temp_68.mat');
% SC_temp_69 = load('projects/axon_simulations/SC_temp_data/SC_temp_69.mat');
% SC_temp_70 = load('projects/axon_simulations/SC_temp_data/SC_temp_70.mat');
% SC_temp_71 = load('projects/axon_simulations/SC_temp_data/SC_temp_71.mat');

% 
% DC_temp_35 = load('projects/axon_simulations/DC_temp_data/DC_temp_35.mat');
% DC_temp_37 = load('projects/axon_simulations/DC_temp_data/DC_temp_37.mat');
% DC_temp_39 = load('projects/axon_simulations/DC_temp_data/DC_temp_39.mat');
% DC_temp_41 = load('projects/axon_simulations/DC_temp_data/DC_temp_41.mat');
% DC_temp_43 = load('projects/axon_simulations/DC_temp_data/DC_temp_43.mat');
% DC_temp_45 = load('projects/axon_simulations/DC_temp_data/DC_temp_45.mat');
% DC_temp_47 = load('projects/axon_simulations/DC_temp_data/DC_temp_47.mat');
% DC_temp_49 = load('projects/axon_simulations/DC_temp_data/DC_temp_49.mat');
% DC_temp_51 = load('projects/axon_simulations/DC_temp_data/DC_temp_51.mat');
% DC_temp_53 = load('projects/axon_simulations/DC_temp_data/DC_temp_53.mat');
% DC_temp_55 = load('projects/axon_simulations/DC_temp_data/DC_temp_55.mat');
% DC_temp_57 = load('projects/axon_simulations/DC_temp_data/DC_temp_57.mat');
% DC_temp_59 = load('projects/axon_simulations/DC_temp_data/DC_temp_59.mat');
% DC_temp_61 = load('projects/axon_simulations/DC_temp_data/DC_temp_61.mat');
% DC_temp_63 = load('projects/axon_simulations/DC_temp_data/DC_temp_63.mat');
% DC_temp_64 = load('projects/axon_simulations/DC_temp_data/DC_temp_64.mat');
% DC_temp_65 = load('projects/axon_simulations/DC_temp_data/DC_temp_65.mat');
% DC_temp_66 = load('projects/axon_simulations/DC_temp_data/DC_temp_66.mat');
% DC_temp_67 = load('projects/axon_simulations/DC_temp_data/DC_temp_67.mat');
% DC_temp_68 = load('projects/axon_simulations/DC_temp_data/DC_temp_68.mat');
% DC_temp_69 = load('projects/axon_simulations/DC_temp_data/DC_temp_69.mat');
% DC_temp_70 = load('projects/axon_simulations/DC_temp_data/DC_temp_70.mat');
% DC_temp_71 = load('projects/axon_simulations/DC_temp_data/DC_temp_71.mat');

% picking time shots
time1 = 5.1; % in ms
time2 = 8; % in ms
time3 = 8.5; % in ms
time4 = 9; % in ms
time5 = 9.5; % in ms
time6 = 10; % in ms
time7 = 10.5; % in ms

% time1 = 1; % in ms
% time2 = 2; % in ms
% time3 = 2.5; % in ms
% time4 = 3; % in ms
% time5 = 3.5; % in ms
% time6 = 4; % in ms
% time7 = 5; % in ms

list_of_times = {time1, time5, time6};

% picking space shots
position1 = 0.5; % in cm
position2 = 1; % in cm
position3 = 1.5; % in cm
position4 = 2; % in cm
position5 = 2.5; % in cm
position6 = 3; % in cm
position7 = 5; % in cm

% position1 = 0.015; % in cm
% position2 = 0.02; % in cm
% position3 = 0.025; % in cm
% position4 = 0.03; % in cm
% position5 = 0.035; % in cm
% position6 = 0.038; % in cm
% position7 = 0.04; % in cm

list_of_positions = {position2, position6, position7};


% picking pause (this controls the speed of the animation, the pause variable 
% is in seconds. It is how many seconds each time (or space) shot will be
% paused at). NOTE: the legends are what slows down the animations, may
% condsider getting rid of or adding them in certain cases.
p = 0.01;

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
% set_of_data10 = {SC_Cohen_Optimized_params, SC_Cohen_Avg_params, SC_Cohen_DC_Optimized_params, DC_Cohen_Optimized_params, DC_Cohen_Avg_params};
% set_of_data11 = {HH_temp_base, HH_temp_8, HH_temp_10, HH_temp_12, HH_temp_14, HH_temp_16, HH_temp_18, HH_temp_20, HH_temp_22, HH_temp_24, HH_temp_26, HH_temp_28, HH_temp_30, HH_temp_31, HH_temp_32, HH_temp_33, HH_temp_34, HH_temp_35};
% set_of_data12 = {SC_temp_35, SC_temp_37, SC_temp_39, SC_temp_41, SC_temp_43, SC_temp_45, SC_temp_47, SC_temp_49, SC_temp_51, SC_temp_53, SC_temp_55, SC_temp_57, SC_temp_59, SC_temp_61, SC_temp_63, SC_temp_65, SC_temp_66, SC_temp_67, SC_temp_68, SC_temp_69, SC_temp_70, SC_temp_71};
% set_of_data13 = {DC_temp_35, DC_temp_37, DC_temp_39, DC_temp_41, DC_temp_43, DC_temp_45, DC_temp_47, DC_temp_49, DC_temp_51, DC_temp_53, DC_temp_55, DC_temp_57, DC_temp_59, DC_temp_61, DC_temp_63, DC_temp_64, DC_temp_65, DC_temp_66, DC_temp_67, DC_temp_68, DC_temp_69, DC_temp_70, DC_temp_71};

% data = SC_Cohen_Optimized_params;
data = HH_1;

plot_animation_voltage_vs_time(data, p);
% plot_animation_voltage_vs_space(data, p);
% plot_animation_probabilities_vs_time(HH_data_Temp_33, p);
% plot_animation_probabilities_vs_space(HH_data_Temp_32, p);
% plot_time_and_space_shots(data, list_of_positions, list_of_times);
% plot_voltage_vs_time_comparison(set_of_data8, p);
% plot_voltage_vs_space_comparison(set_of_data12, p);
% plot_Vm_and_Vm_minus_Vmy_vs_space(data, p)



%%%%%%%%%%%%%%%%%%%%%
% PLOTTER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%

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
    
    xmin = 0;
    xmax = T;
    ymin = -70;
    ymax = 50;
    
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
        legendStrings1{end+1} = sprintf('$V_m$ at t = %g ms', list_of_times{i});
    end
    legend(legendStrings1, 'Interpreter','latex')
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel("Length of the axon in cm.")
    
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
        legendStrings2{end+1} = sprintf('$V_m$ at x = %g cm', list_of_positions{i});
    end
    legend(legendStrings2, 'Interpreter', 'latex')
    
    xmin = 0;
    xmax = T;
    ymin = -80;
    ymax = 50;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel("Time in milliseconds.")
    
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
    legend(legendStrings3, 'Interpreter','latex')
    ylabel("Probabilities of ion channels opening/closing.")
    xlabel("Time in milliseconds.")
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
                plot(t, data_set{j}.(voltages{k})(:,i), 'Color', [1, 0, 0] * j/length(data_set));
                hold on
            end

            % legend('SiGe Tube params', 'Tube+Paralyne params', 'Location', 'northeast');
            legend('SC model: $R_i = 0.144 k\Omega cm$', 'DC model: $R_i = 0.712 k\Omega cm$', 'Location', 'northeast', 'Interpreter', 'latex')
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
                plot(t, data_set{j}.(voltages{k})(i,:), 'Color', [1, 0, 0] * j/length(data_set));
                hold on
            end

            % Add the legend (NOTE: the legend is what is slowing down the animation)
            legend('DC: Myelinated', 'DC: SiGe Tube params', 'DC: Tube+Paralyne params', 'Location', 'northeast');
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
        title('4) DC: cell 6 params Spatial Profile', 'Interpreter', 'latex');

        % Add the legend (NOTE: the legend is what is slowing down the animation)
        legend('$V_m$', '$V_m - V_{my}$', 'Location', 'northeast', 'Interpreter', 'latex');
        text(xmin+0.01, ymax-10, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
    
        % Add a pause to create animation effect
        pause(p);

    end
end


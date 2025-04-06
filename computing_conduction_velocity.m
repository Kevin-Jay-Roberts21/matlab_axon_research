% A MATLAB script to compute the conduction velocity of the voltage
% simulations. More accurate results for longer simulations.
% Kevin Roberts
% February 2025

clear all
close all
clc

% Collecting temp data
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
% SC_temp_20 = load('projects/axon_simulations/SC_temp_data/SC_temp_base.mat');
% SC_temp_22 = load('projects/axon_simulations/SC_temp_data/SC_temp_22.mat');
% SC_temp_24 = load('projects/axon_simulations/SC_temp_data/SC_temp_24.mat');
% SC_temp_26 = load('projects/axon_simulations/SC_temp_data/SC_temp_26.mat');
% SC_temp_28 = load('projects/axon_simulations/SC_temp_data/SC_temp_28.mat');
% SC_temp_30 = load('projects/axon_simulations/SC_temp_data/SC_temp_30.mat');
% SC_temp_32 = load('projects/axon_simulations/SC_temp_data/SC_temp_32.mat');
% SC_temp_34 = load('projects/axon_simulations/SC_temp_data/SC_temp_34.mat');
% SC_temp_36 = load('projects/axon_simulations/SC_temp_data/SC_temp_36.mat');
% SC_temp_38 = load('projects/axon_simulations/SC_temp_data/SC_temp_38.mat');
% SC_temp_40 = load('projects/axon_simulations/SC_temp_data/SC_temp_40.mat');
% SC_temp_42 = load('projects/axon_simulations/SC_temp_data/SC_temp_42.mat');
% SC_temp_44 = load('projects/axon_simulations/SC_temp_data/SC_temp_44.mat');
% SC_temp_46 = load('projects/axon_simulations/SC_temp_data/SC_temp_46.mat');
% SC_temp_48 = load('projects/axon_simulations/SC_temp_data/SC_temp_48.mat');
% SC_temp_50 = load('projects/axon_simulations/SC_temp_data/SC_temp_50.mat');
% SC_temp_52 = load('projects/axon_simulations/SC_temp_data/SC_temp_52.mat');
% SC_temp_54 = load('projects/axon_simulations/SC_temp_data/SC_temp_54.mat');
% SC_temp_55 = load('projects/axon_simulations/SC_temp_data/SC_temp_55.mat');
% 
% DC_temp_20 = load('projects/axon_simulations/DC_temp_data/DC_temp_base.mat');
% DC_temp_21 = load('projects/axon_simulations/DC_temp_data/DC_temp_21.mat');
% DC_temp_22 = load('projects/axon_simulations/DC_temp_data/DC_temp_22.mat');
% DC_temp_23 = load('projects/axon_simulations/DC_temp_data/DC_temp_23.mat');
% DC_temp_24 = load('projects/axon_simulations/DC_temp_data/DC_temp_24.mat');


% SC_data_from_saltcond2023_code = load('SC_data_from_saltcond2023_code.mat');
% DC_data_from_saltcond2023_code = load('DC_data_from_saltcond2023_code.mat');

SC_Huang_Myelinated = load('projects/axon_simulations/Huang_simulations/SC_Huang_Myelinated.mat');
SC_Huang_Tube_params = load('projects/axon_simulations/Huang_simulations/SC_Huang_Tube_params.mat');
SC_Huang_TubeParalyne_params = load('projects/axon_simulations/Huang_simulations/SC_Huang_TubeParalyne_params.mat');
% 
% SC_Cohen_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_cell6_params.mat');
% SC_Cohen_avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_avg_params.mat');
% SC_Cohen_DC_avg_params = load('projects/axon_simulations/Cohen_param_simulations/SC_Cohen_DC_avg_params.mat');
% DC_Cohen_cell6_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_cell6_params.mat');
% DC_Cohen_avg_params = load('projects/axon_simulations/Cohen_param_simulations/DC_Cohen_avg_params.mat');

% HH_1 = load('HH_1.mat');

% NOTE: Time and space choices will vary depending on HH, SC, and DC models

% creating a list of spaces (x1 and x2) to compute the cv, we will then
% take an average of these different spaces to get a more accurate cv
% NOTE: these may differ for different length axons
HH_space_vec = [4 5; 5 6; 6 7; 7 8;]; % all in (cm)
HH_space_vec1 = [2 2.5; 2.5 3; 3 3.5; 3.5 4];
SC_and_DC_space_vec = [0.06 0.08; 0.08 0.1; 0.1 0.12; 0.12 0.14;]; % all in (cm)

% defining multiple datasets
% HH_data_set = {HH_data_Temp_base, HH_data_Temp_7, HH_data_Temp_8, HH_data_Temp_9, HH_data_Temp_10, HH_data_Temp_11, HH_data_Temp_12, HH_data_Temp_13, HH_data_Temp_14, HH_data_Temp_15, HH_data_Temp_16, HH_data_Temp_17, HH_data_Temp_18, HH_data_Temp_19, HH_data_Temp_20, HH_data_Temp_21, HH_data_Temp_22, HH_data_Temp_23, HH_data_Temp_24, HH_data_Temp_25, HH_data_Temp_26, HH_data_Temp_27, HH_data_Temp_28, HH_data_Temp_29, HH_data_Temp_30, HH_data_Temp_31, HH_data_Temp_32};
% SC_data_set = {SC_data_Temp_20, SC_data_Temp_22, SC_data_Temp_24, SC_data_Temp_26, SC_data_Temp_28, SC_data_Temp_30, SC_data_Temp_32, SC_data_Temp_34, SC_data_Temp_36, SC_data_Temp_38, SC_data_Temp_40, SC_data_Temp_42, SC_data_Temp_44, SC_data_Temp_46, SC_data_Temp_48, SC_data_Temp_50};
% DC_data_set = {DC_data_Temp_20, DC_data_Temp_22, DC_data_Temp_24, DC_data_Temp_26, DC_data_Temp_28, DC_data_Temp_30, DC_data_Temp_32, DC_data_Temp_34, DC_data_Temp_36, DC_data_Temp_38, DC_data_Temp_40, DC_data_Temp_42, DC_data_Temp_44, DC_data_Temp_46, DC_data_Temp_48, DC_data_Temp_50};
% data_set1 = {SC_data_from_saltcond2023_code, DC_data_from_saltcond2023_code};
data_set2 = {SC_Huang_Myelinated, SC_Huang_Tube_params, SC_Huang_TubeParalyne_params};
% data_set3 = {SC_Cohen_cell6_params, SC_Cohen_avg_params, SC_Cohen_DC_avg_params, DC_Cohen_cell6_params, DC_Cohen_avg_params};
% data_set4 = {};
% data_set5 = {HH_temp_base, HH_temp_8, HH_temp_10, HH_temp_12, HH_temp_14, HH_temp_16, HH_temp_18, HH_temp_20, HH_temp_22, HH_temp_24, HH_temp_26, HH_temp_28, HH_temp_30, HH_temp_31, HH_temp_32};
% data_set12 = {SC_temp_20, SC_temp_22, SC_temp_24, SC_temp_26, SC_temp_28, SC_temp_30, SC_temp_32, SC_temp_34, SC_temp_36, SC_temp_38, SC_temp_40, SC_temp_42, SC_temp_44, SC_temp_46, SC_temp_48, SC_temp_50, SC_temp_52, SC_temp_54, SC_temp_55};
% data_set7 = {DC_temp_20, DC_temp_21, DC_temp_22, DC_temp_23};

% picking the dataset to compute cv
data = data_set2;

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




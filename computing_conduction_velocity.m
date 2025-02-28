% A MATLAB script to compute the conduction velocity of the voltage
% simulations. More accurate results for longer simulations.
% Kevin Roberts
% February 2025

clear all
close all
clc

% Collecting data
% getting the saved data
% HH_data = load('HH_data.mat');
HH_data1 = load('HH_1.mat');
% HH_data2 = load('HH_data2.mat');
% SC_data = load('SC_data.mat');
% SC_data1 = load('SC_data1.mat');
% SC_data2 = load('SC_data2.mat');

% Collecting temp data
% HH_data_Temp_default = load('axon_simulations/HH_temp_data/HH_data_Temp_6.3.mat');
% HH_data_Temp_7 = load('axon_simulations/HH_temp_data/HH_data_Temp_7.mat');
% HH_data_Temp_8 = load('axon_simulations/HH_temp_data/HH_data_Temp_8.mat');
% HH_data_Temp_9 = load('axon_simulations/HH_temp_data/HH_data_Temp_9.mat');
% HH_data_Temp_10 = load('axon_simulations/HH_temp_data/HH_data_Temp_10.mat');
% HH_data_Temp_15 = load('axon_simulations/HH_temp_data/HH_data_Temp_15.mat');
% HH_data_Temp_20 = load('axon_simulations/HH_temp_data/HH_data_Temp_20.mat');
% HH_data_Temp_25 = load('axon_simulations/HH_temp_data/HH_data_Temp_25.mat');
% HH_data_Temp_30 = load('axon_simulations/HH_temp_data/HH_data_Temp_30.mat');

% To get the conduction velocity, we pick a space vector vec_x1 at some time 
% t1, and calculate the position at where the max of this vector is, call 
% it x1. Then we pick another space vector vec_x2 at some time t2 where 
% t2 > t1, and calculate the position where the max of vec_x2 is, call it
% x2. Then the conduction velocity is calculated as cv = (x2 - x1)/(t2 - t1)

% NOTE: Time and space choices will vary depending on HH, SC, and DC models

% picking the times:
t1 = 10; % (in ms)
t2 = 10.3; % (in ms)

% picking the dataset to compute cv
data = HH_data1;

% computing conduction velocity
cv1 = calculate_cv(t1, t2, data) % (in cm/ms)

% converting cv to m/s
cv2 = cv1 * 10


% conduction velocity function
function cv = calculate_cv(t1, t2, data)
    
    % need to find the index of t1 and t2
    index_t1 = t1/data.dt;
    index_t2 = t2/data.dt;
    
    % identifying the space vectors at index_t1 and index_t2
    vec_x1 = data.Vm_all(index_t1,:);
    vec_x2 = data.Vm_all(index_t2,:);

    % computing the positions of where the max voltage is at t1 and t2:
    [Vm_at_x1, index_x1] = max(vec_x1); % returns max voltage at t1 and index of where the max is
    [Vm_at_x2, index_x2] = max(vec_x2); % returns max voltage at t2 and index of where the max is
    
    % Note that x1 and x2 are just the indices. We need to find their actual
    % spatial position in cm. This is done by identifiying the mesh from the
    % data
    x1 = index_x1*data.dx; % (in cm)
    x2 = index_x2*data.dx; % (in cm)
    
    % finally, compute the conduction velocity
    cv = (x2 - x1)/(t2 - t1);

end




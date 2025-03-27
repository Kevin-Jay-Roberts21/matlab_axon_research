% A script to used to compute and plot the time it takes to generate an
% action potential given a stimulus value (Only for HH equation)
% Kevin Roberts
% March 2025

clear all
close all
clc

% uploading the data 
HH_S_100 = load('projects/axon_simulations/HH stim threshold/HH_S_100.mat');
HH_S_75 = load('projects/axon_simulations/HH stim threshold/HH_S_75.mat');
HH_S_50 = load('projects/axon_simulations/HH stim threshold/HH_S_50.mat');
HH_S_25 = load('projects/axon_simulations/HH stim threshold/HH_S_25.mat');
HH_S_15 = load('projects/axon_simulations/HH stim threshold/HH_S_15.mat');
HH_S_14 = load('projects/axon_simulations/HH stim threshold/HH_S_14.mat');
HH_S_13 = load('projects/axon_simulations/HH stim threshold/HH_S_13.mat');
HH_S_12 = load('projects/axon_simulations/HH stim threshold/HH_S_12.mat');
HH_S_11_5 = load('projects/axon_simulations/HH stim threshold/HH_S_11.5.mat');
HH_S_11_41 = load('projects/axon_simulations/HH stim threshold/HH_S_11.41.mat');
HH_S_11_406 = load('projects/axon_simulations/HH stim threshold/HH_S_11.406.mat');
HH_S_11_4055 = load('projects/axon_simulations/HH stim threshold/HH_S_11.4055.mat');
HH_S_11_40546 = load('projects/axon_simulations/HH stim threshold/HH_S_11.40546.mat');
HH_S_11_405456 = load('projects/axon_simulations/HH stim threshold/HH_S_11.405456.mat');
HH_S_11_4054556 = load('projects/axon_simulations/HH stim threshold/HH_S_11.4054556.mat');
HH_S_11_40545555 = load('projects/axon_simulations/HH stim threshold/HH_S_11.40545555.mat');
HH_S_11_405455545 = load('projects/axon_simulations/HH stim threshold/HH_S_11.405455545.mat');
HH_S_11_4054555448 = load('projects/axon_simulations/HH stim threshold/HH_S_11.4054555448.mat');

list_of_data = {HH_S_100, HH_S_50, HH_S_12, HH_S_11_5, HH_S_11_41, HH_S_11_406, HH_S_11_4055, HH_S_11_40546, HH_S_11_405456, HH_S_11_4054556, HH_S_11_40545555, HH_S_11_405455545, HH_S_11_4054555448};

% creating the stimulus grid
stim_grid = [100, 50, 12, 11.5, 11.41, 11.406, 11.4055, 11.40546, 11.405456, 11.4054556, 11.40545555, 11.405455545, 11.4054555448];

% creating a list of times between the stimulus and the peak of the action
% potential for each of the data sets given
action_potential_times = zeros(1, length(stim_grid));
time_after_stimulus_interval = HH_S_100.S_T1; % (ms)
space_after_stimulus_interval = HH_S_100.S_P0; % (cm)
dx = HH_S_100.dx;
dt = HH_S_100.dt;

vec = HH_S_100.Vm_all(:,1);

for i = 1:length(list_of_data)
    
    % getting the time vector right after stimulus is added
    time_vec = list_of_data{i}.Vm_all(:,space_after_stimulus_interval/dx);
    
    % getting the max of the time vector (which is the time when the peak 
    % of action potential occurs) after the stimulus is added
    [voltage_of_AP_peak, index_of_AP_peak] = max(time_vec);
    
    % calculating the time from stimulus to AP peak
    time_diff = index_of_AP_peak*dt - time_after_stimulus_interval;
    
    % Adding it to the vector
    action_potential_times(i) = time_diff;

end 

% Example data
x = 0:0.1:10;
y = sin(x);

% Create the plot
figure;
plot(x, y);

% Fix y-axis limits
ylim([-1.5, 1.5]);

% Set y-tick positions (e.g., at -1, 0, 1)
yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);

% Set custom tick labels
yticklabels({'0', '0', '0', '3.5', '4.5', '2', '1'});

% Optionally, add labels to the plot
xlabel('X axis');
ylabel('Y axis');
title('Custom Y-axis Ticks');
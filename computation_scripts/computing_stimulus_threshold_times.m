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
list_of_data1 = {HH_S_12, HH_S_11_5, HH_S_11_41, HH_S_11_406, HH_S_11_4055, HH_S_11_40546, HH_S_11_405456, HH_S_11_4054556, HH_S_11_40545555, HH_S_11_405455545, HH_S_11_4054555448};
list_of_data2 = {HH_S_12, HH_S_11_406, HH_S_11_4054556};

% creating the stimulus grid
stim_grid = [100, 50, 12, 11.5, 11.41, 11.406, 11.4055, 11.40546, 11.405456, 11.4054556, 11.40545555, 11.405455545, 11.4054555448];
stim_grid1 = [12, 11.5, 11.41, 11.406, 11.4055, 11.40546, 11.405456, 11.4054556, 11.40545555, 11.405455545, 11.4054555448];
stim_grid2 = [12, 11.406, 11.4054556];

% creating a list of times between the stimulus and the peak of the action
% potential for each of the data sets given
action_potential_times = zeros(1, length(stim_grid1));
time_after_stimulus_interval = HH_S_100.S_T1; % (ms)
space_after_stimulus_interval = HH_S_100.S_P0; % (cm)
dx = HH_S_100.dx;
dt = HH_S_100.dt;

vec = HH_S_100.Vm_all(:,1);

for i = 1:length(list_of_data1)

    % getting the time vector right after stimulus is added
    time_vec = list_of_data1{i}.Vm_all(:,space_after_stimulus_interval/dx);

    % getting the max of the time vector (which is the time when the peak 
    % of action potential occurs) after the stimulus is added
    [voltage_of_AP_peak, index_of_AP_peak] = max(time_vec);

    % calculating the time from stimulus to AP peak
    time_diff = index_of_AP_peak*dt - time_after_stimulus_interval;

    % Adding it to the vector
    action_potential_times(i) = time_diff;

end 

% Compute difference from threshold
threshold = min(stim_grid1);
delta_stim = stim_grid1 - threshold;

figure;
semilogx(delta_stim, action_potential_times, 'o-', 'LineWidth', 2);
xlabel('Stimulus above threshold (log scale)');
ylabel('Action Potential Time (ms)');
grid on;












% % Flip x-axis
% set(gca, 'XDir','reverse');
% 
% % Example data
% action_potential_times = fliplr(action_potential_times);
% y = [11.4054555448, 18.78833425, 26.17121295, 33.55409166, 40.93697036, 48.31984907, 55.70272777, 63.08560648, 70.46848518, 77.85136389, 85.23424259, 92.6171213, 100];
% 
% % Manually spaced stimulus values (not to scale)
% stim_values = [12, 11.45, 11.4054556];  % Slide 11.406 more left visually
% ap_times = fliplr(action_potential_times);  % Make sure this matches the 3 values
% 
% % Plot
% figure;
% plot(stim_values, ap_times, '-o', ...
%     'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
% 
% % Ensure x-axis increases left to right (12 on the right)
% set(gca, 'XDir','normal');
% 
% % Label axes
% xlabel('Stimulus value (k$\Omega$ cm$^2$)', 'Interpreter','latex');
% ylabel('Time to AP peak (ms)', 'Interpreter','latex');
% 
% % Define the new 'stim. threshold' tick further left
% x_tick_stim_thresh = 11.35;  % Moved further left
% xticks_all = sort([x_tick_stim_thresh, stim_values]);
% 
% % Define custom labels (must match number of ticks)
% xticklabels_all = {'stim. threshold', '11.4054556', '11.406', '12'};
% 
% % Apply ticks and labels
% xticks(xticks_all);
% xticklabels(xticklabels_all);
% 
% % Slant x-axis labels to avoid overlap
% xtickangle(45);  % You can adjust the angle as needed
% 
% % Adjust limits and grid
% xlim([11.3, 12.1]);  % Leave room for the new label
% ylim([0, max(ap_times)*1.1]);
% grid on;
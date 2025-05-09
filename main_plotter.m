% Script created to use plotter functions to plot data collected from
% running simulations
% Kevin Roberts
% May 2025

clear all
close all
clc

% There are 9 total different types of data that we would like plot. We
% will categorize them as follows:
% data_type1: "Vm", "Vmy", "Vm_minus_Vmy", "n", "m" and "h"
% data_type2: "nmh", "Vm_and_Vm_minus_Vmy", and
% "Vm_and_Vmy_and_Vm_minus_Vmy"
% data_type3: everything in data_type1 AND data_type2

% We note that in data_type2, multiple lines are plotted, whereas datatype1
% we only one line is plotted. 


% Loading in the data to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: some data may be saved in different paths
hh_simulation = load(pwd + "/hh_simulation.mat");
sc_simulation_v3 = load(pwd + "/sc_simulation_v3.mat");
dc_simulation_v3 = load(pwd + "/dc_simulation_v3.mat");

%%
% The following variables are subject to change due to users preference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining data_set
data_set = {sc_simulation_v3, dc_simulation_v3};
data_set_names = {'sc_simulation_v3', 'dc_simulation_v3'};

% Defining data
data = sc_simulation_v3;

% Defining p (controls animation speed, smaller p implies faster animation)
p = 0.01;

% Defining the set of space shots
% position1 = 1; % in cm
% position2 = 3; % in cm
% position3 = 5; % in cm
space_shot1 = 0.0005; % in cm
space_shot2 = 0.06; % in cm
space_shot3 = 0.16; % in cm
space_shots = {space_shot1, space_shot2, space_shot3};

% Defining the set of time shots
% time_shot1 = 5.1; % in ms
% time_shot2 = 9.5; % in ms
% time_shot3 = 10; % in ms
time_shot1 = 2; % in ms
time_shot2 = 3.3; % in ms
time_shot3 = 4; % in ms
time_shots = {time_shot1, time_shot2, time_shot3};

% Defining a single space shot
space_shot = 0.08; % in cm

% Defining a single time shot
time_shot = 1.1; % in ms

% Defining axon segment
position1 = 0.0800; % in cm
position2 = 0.0828; % in cm
position3 = 0.0885; % in cm
axon_segment = {position1, position2, position3};

% Defining multiple axon segments
axon_segment1 = [0.0400 0.0485]; % interval is in cm
axon_segment2 = [0.0800 0.0885]; % interval is in cm
axon_segment3 = [0.1200 0.1285];
axon_segments = {axon_segment1, axon_segment2, axon_segment3};


% Choosing a data_type to plot
data_type1 = "Vmy"; % can choose any other variable in data_type1 set
data_type2 = "Vm_and_Vm_minus_Vmy"; % can choose any other variable in data_type2 set
data_type3 = "Vm_and_Vmy_and_Vm_minus_Vmy"; % can choose any other variable in data_type3 set



%%%%%%%%%%%%%%%%
% All Plotters %
%%%%%%%%%%%%%%%%

% Animations
%%%%%%%%%%%%
%%
plotter_functions.plot_animation_data_vs_time(data, data_type3, p);

%%
plotter_functions.plot_animation_data_vs_space(data, data_type3, p);

%%
plotter_functions.plot_animation_comparison_data_vs_time(data_set, data_set_names, data_type1, p);

%%
plotter_functions.plot_animation_comparison_data_vs_space(data_set, data_set_names, data_type1, p);


% Time and Space Shots
%%%%%%%%%%%%%%%%%%%%%%
%%
plotter_functions.plot_data_vs_time_at_space_shots(data, data_type1, space_shots);

%%
plotter_functions.plot_data_vs_space_at_time_shots(data, data_type1, time_shots);

%%
plotter_functions.plot_data_vs_time_at_one_space_shot(data, data_type3, space_shot);

%%
plotter_functions.plot_data_vs_space_at_one_time_shot(data, data_type3, time_shot);

%%
plotter_functions.plot_data_vs_space_axon_segment(data, data_type1, axon_segment);

%%
plotter_functions.plot_data_vs_space_axon_segments(data, data_type1, axon_segments);


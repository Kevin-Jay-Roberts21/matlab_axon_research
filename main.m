% main matlab script used to run simulations, collect data and plot results
% Kevin Roberts
% April 2025

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running an HH Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting the mesh and material parameters
hh_mesh = hh_mesh_parameters();
hh_material = hh_material_parameters();

% Running the HH simulation
hh_simultion = hh_function_v1(hh_mesh, hh_material);

% % Saving HH data
% save('hh_simulation.mat', 'hh_simultion');


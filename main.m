% main matlab script used to run simulations, collect data and plot results
% Kevin Roberts
% April 2025

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running an HH Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Getting the mesh and material parameters
% hh_mesh = hh_mesh_parameters();
% hh_material = hh_material_parameters();
% 
% % Running the HH simulation
% hh_simultion = hh_function_v1(hh_mesh, hh_material);
% 
% % % Saving HH data
% % save('hh_simulation.mat', 'hh_simultion');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running a SC Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting the mesh and material parameters
sc_mesh = sc_dc_mesh_parameters();
sc_material = set1_parameters();

% Running the SC simulation
sc_simulation = sc_function_v2(sc_mesh, sc_material);

% % Saving SC data
% save('sc_simulation.mat', 'sc_simultion');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running a DC Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Getting the mesh and material parameters
% dc_mesh = sc_dc_mesh_parameters();
% dc_material = set1_parameters();
% 
% % Running the SC simulation
% dc_simulation = dc_function_v3(dc_mesh, dc_material);
% 
% % % Saving SC data
% % save('sc_simulation.mat', 'sc_simultion');




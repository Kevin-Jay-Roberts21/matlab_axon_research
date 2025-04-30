% main matlab script used to run simulations, collect data and plot results
% Kevin Roberts
% April 2025

clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running an HH Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the code below to run an hh simulation (data will automatically be saved)

% Getting the mesh and material parameters
hh_mesh = mesh_parameter_sets.hh_mesh_parameters();
hh_material = material_parameter_sets.hh_material_parameters();

% Running the HH simulation
hh_simultion = numerical_scheme_functions.hh_function_v1(hh_mesh, hh_material);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running a SC Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the code below to run a sc simulation (data will automatically be saved)

% Getting the mesh and material parameters
sc_mesh = mesh_parameter_sets.sc_dc_mesh_parameters();
sc_material = material_parameter_sets.set1_parameters();

% Running the SC simulation
sc_simulation = numerical_scheme_functions.sc_function_v1(sc_mesh, sc_material);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running a DC Simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the code below to run a dc simulation (data will automatically be saved)


% Getting the mesh and material parameters
dc_mesh = mesh_parameter_sets.sc_dc_mesh_parameters();
dc_material = material_parameter_sets.set1_parameters();

% Running the SC simulation
dc_simulation = numerical_scheme_functions.dc_function_v3(dc_mesh, dc_material);




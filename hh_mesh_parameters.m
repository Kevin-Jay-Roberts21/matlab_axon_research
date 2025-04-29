% HH defualt mesh parameters
% Kevin Roberts
% April 2025


function hh_mesh_params = hh_mesh_parameters()
    % Defining the mesh parameters
    hh_mesh_params.a = 0.025; % (cm) axon radius 
    hh_mesh_params.L = 5; % (cm) axon length 
    hh_mesh_params.dx = 0.01; % (cm) space step 
    hh_mesh_params.T = 30; % (ms) total time 
    hh_mesh_params.dt = 0.01; % (ms) time step 
    hh_mesh_params.m = L/dx + 1; % (#) total number of space steps 
    hh_mesh_params.n = T/dt + 1; % (#) total number of time steps 
end
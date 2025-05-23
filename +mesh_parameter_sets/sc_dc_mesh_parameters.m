% SC and DC mesh parameters
% Kevin Roberts
% April 2025

function sc_dc_mesh_params = sc_dc_mesh_parameters()

    % Defining the mesh parameters
    sc_dc_mesh_params.dx = 0.00005; % (cm) space step
    sc_dc_mesh_params.dt = 0.01; % (ms) time step 
    sc_dc_mesh_params.L_my = 0.0075; % (cm) internodal length
    sc_dc_mesh_params.L_n = 0.0005; % (cm) nodal length
    sc_dc_mesh_params.L_pn = 2.3*10^(-4); % (cm) paranodal length
    sc_dc_mesh_params.L_s = sc_dc_mesh_params.L_n + sc_dc_mesh_params.L_my; % (cm) length of an axon segment
    sc_dc_mesh_params.n_s = 20; % (#) number of axon segments
    sc_dc_mesh_params.L = sc_dc_mesh_params.n_s*sc_dc_mesh_params.L_s; % (cm) total length of axon
    sc_dc_mesh_params.T = 30; % (ms) the total time of the experiment
    sc_dc_mesh_params.N_n = round(sc_dc_mesh_params.L_n/sc_dc_mesh_params.dx); % (#) number of space steps in a nodal region
    sc_dc_mesh_params.N_my = round(sc_dc_mesh_params.L_my/sc_dc_mesh_params.dx); % (#) number of space steps in an internodal region
    sc_dc_mesh_params.N_s = sc_dc_mesh_params.N_n + sc_dc_mesh_params.N_my; % (#) number of space steps in an entire axon segement
    sc_dc_mesh_params.m = sc_dc_mesh_params.N_s*sc_dc_mesh_params.n_s + 1; % total number of space steps
    sc_dc_mesh_params.n = sc_dc_mesh_params.T/sc_dc_mesh_params.dt + 1; % n is the number of time steps
    
    % NOTE: d_pa, d_pn, and L_pn are not used in the SC model

end
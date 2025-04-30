% HH material parameters
% Kevin Roberts
% April 2025

function hh_material_params = hh_material_parameters()
    
    % Defining all of the material and intrinsic parameters
    hh_material_params.a = 0.025; % (cm) axon radius 
    hh_material_params.C_m = 1; % (micro-farads/cm^2) membrane capacitance 
    hh_material_params.R_i = 0.030; % (kilo-ohms * cm) specific intracellular resistivity 
    hh_material_params.G_K = 36; % (mS/cm^2) maximum potassium conductance
    hh_material_params.G_Na = 120; % (mS/cm^2) maximum sodium conductance
    hh_material_params.G_L = 0.3; % (mS/cm^2) maximum leak conductance
    hh_material_params.E_K = -77; % (mV) Nernst potential for potassium ions 
    hh_material_params.E_Na = 50; % (mV) Nernst potential for sodium ions 
    hh_material_params.E_L = -54.4; % (mV) Nernst potential for leak channels 
    
    % Defining the temperature parameters
    hh_material_params.T_base = 6.3; % (C) base temperature
    hh_material_params.T_actual = 6.3; % (C) temperature of the squid axon fiber
    hh_material_params.Q_10 = 3; % (#) temperature coefficient
    
    % Defining the stimulus value, and time and space interval it is applied to
    hh_material_params.S_v = 20; % (mS/cm^2) % stimulus value
    hh_material_params.S_T0 = 5; % (ms) start time of when stimulus is added 
    hh_material_params.S_T1 = 5.1; % (ms) end time of when stimulus is added 
    hh_material_params.S_P0 = 1; % (cm) start position of adding the stimulus 
    hh_material_params.S_P1 = 1.1; % (cm) end position of adding the stimulus 
    
    % Defining the initial Vm, n, m, h values (uniformly distributed along the axon)
    hh_material_params.N_0 = 0.3177; % (#) initial value of potassium activation gating variable
    hh_material_params.M_0 = 0.0529; % (#) initial value of sodium activation gating variable
    hh_material_params.H_0 = 0.5961; % (#) initial value of sodium inactivation gating
    hh_material_params.V_m0 = -64.9997; % (mV) initial value of membrane voltage

end


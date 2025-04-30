% cell 5 material parameters
% Kevin Roberts
% April 2025

%%%%%%%%%%%%%%%%%%%%%%%%%
% cell 5 Parameter Set %
%%%%%%%%%%%%%%%%%%%%%%%%%
function cell5_params = cell5_parameters()

    % Defining the material properties
    cell5_params.a = 0.55*10^(-4); % (cm) radius in nodal region
    cell5_params.a_my = cell5_params.a/0.698; % (cm) radius in myelinated region
    cell5_params.R_i = 0.109; % (kilo-ohms*cm) intracellular resistivity
    cell5_params.R_m = 27.7; % (kilo-ohms*cm^2) specific membrane resistance
    cell5_params.C_m = 1.01; % (micro-farads/cm^2) specific membrane capacitance
    cell5_params.r_pa = 105*10^6; % (kilo-ohms/cm) periaxonal resistivity per unit length
    cell5_params.r_pn = 877*10^6; % (kilo-ohms/cm) paranodal resitance per unit length (used in BC since r_bar_pn = r_pn * L_pn) 
    cell5_params.R_my = 75.2; % (kilo-ohms*cm^2) specfic myelin resistance
    cell5_params.C_my = 0.0678; % (micro-fards/cm^2) specific myelin capacitance
    cell5_params.R_pa = cell5_params.r_pa*pi*cell5_params.d_pa*(2*cell5_params.a + cell5_params.d_pa); % (kilo-ohms*cm) resistivity of the periaxonal space (computed)
    cell5_params.R_pn = cell5_params.r_pn*pi*cell5_params.d_pn*(2*cell5_params.a + cell5_params.d_pn); % (kilo-ohms*cm) resistivity of the paranodal space (computed)
    cell5_params.G_K = 80; % (mS/cm^2) max specific potassium conductance
    cell5_params.G_Na = 3000; % (mS/cm^2) max specific sodium conductance 
    cell5_params.G_L = 80; % (mS/cm^2) specific leak conductance
    cell5_params.E_K = -82; % (mV) Nernst potential for potassium ions
    cell5_params.E_Na = 45; % (mV) Nernst potential for sodium ions
    cell5_params.E_L = -59.4; % (mV) Nernst potential for leak channels
    cell5_params.E_rest = -59.4; % (mV) effective resting nernst potential
    
    % Defining the stimulus value, and time and space interval it is applied to
    cell5_params.S_v = 1000; % (mS/cm^2) stimulus value
    cell5_params.S_T0 = 1; % (ms) start time of when stimulus is added 
    cell5_params.S_T1 = 1.1; % (ms) end time of when stimulus is added 
    cell5_params.S_P0 = 0.0001; % (cm) start position of adding the stimulus 
    cell5_params.S_P1 = 0.0004; % (cm) end position of adding the stimulus 
    
    % Defining the temperature parameters
    cell5_params.T_base = 20; % (C) base temperature
    cell5_params.T_actual = 20; % (C) temperature of the squid axon fiber
    cell5_params.Q_10_Na = 2.2; % (#) temperature coefficient for Na current
    cell5_params.Q_10_K = 3; % (#) temperature coefficient for K current 
    
    % Defining the initial Vm, Vmy, n, m, h values (uniformly distributed along the axon)
    cell5_params.V_m0 = -58.1124; % (mV) initial condition for membrane potential 
    cell5_params.V_my0 = 1.2727; % (mV) initial condition for axon potential in periaxonal space
    cell5_params.N_0 = 0.4264; % (#) initial condition for gating variable n
    cell5_params.M_0 = 0.1148; % (#) initial condition for gating variable m
    cell5_params.H_0 = 0.3548; % (#) initial condition for gating variable h

end
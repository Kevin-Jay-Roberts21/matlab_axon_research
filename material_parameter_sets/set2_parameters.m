% Set (2) material parameters
% Kevin Roberts
% April 2025

%%%%%%%%%%%%%%%%%%%%%%%%%
% Set (2) Parameter Set %
%%%%%%%%%%%%%%%%%%%%%%%%%
function set2_params = set2_parameters()

    % Defining the material properties
    set2_params.a = 0.55*10^(-4); % (cm) radius in nodal region
    set2_params.a_my = set2_params.a/0.698; % (cm) radius in myelinated region
    set2_params.R_i = 0.155; % (kilo-ohms*cm) intracellular resistivity
    set2_params.R_m = 24.6; % (kilo-ohms*cm^2) specific membrane resistance
    set2_params.C_m = 1.15; % (micro-farads/cm^2) specific membrane capacitance
    set2_params.r_pa = 125*10^6; % (kilo-ohms/cm) periaxonal resistivity per unit length
    set2_params.r_pn = 2450*10^6; % (kilo-ohms/cm) paranodal resitance per unit length (used in BC since r_bar_pn = r_pn * L_pn) 
    set2_params.R_my = 240; % (kilo-ohms*cm^2) specfic myelin resistance
    set2_params.C_my = 0.0379; % (micro-fards/cm^2) specific myelin capacitance
    set2_params.R_pa = set2_params.r_pa*pi*set2_params.d_pa*(2*set2_params.a + set2_params.d_pa); % (kilo-ohms*cm) resistivity of the periaxonal space (computed)
    set2_params.R_pn = set2_params.r_pn*pi*set2_params.d_pn*(2*set2_params.a + set2_params.d_pn); % (kilo-ohms*cm) resistivity of the paranodal space (computed)
    set2_params.G_K = 80; % (mS/cm^2) max specific potassium conductance
    set2_params.G_Na = 3000; % (mS/cm^2) max specific sodium conductance 
    set2_params.G_L = 80; % (mS/cm^2) specific leak conductance
    set2_params.E_K = -82; % (mV) Nernst potential for potassium ions
    set2_params.E_Na = 45; % (mV) Nernst potential for sodium ions
    set2_params.E_L = -59.4; % (mV) Nernst potential for leak channels
    set2_params.E_rest = -59.4; % (mV) effective resting nernst potential
    
    % Defining the stimulus value, and time and space interval it is applied to
    set2_params.S_v = 1000; % (mS/cm^2) stimulus value
    set2_params.S_T0 = 1; % (ms) start time of when stimulus is added 
    set2_params.S_T1 = 1.1; % (ms) end time of when stimulus is added 
    set2_params.S_P0 = 0.0001; % (cm) start position of adding the stimulus 
    set2_params.S_P1 = 0.0004; % (cm) end position of adding the stimulus 
    
    % Defining the temperature parameters
    set2_params.T_base = 20; % (C) base temperature
    set2_params.T_actual = 20; % (C) temperature of the squid axon fiber
    set2_params.Q_10_Na = 2.2; % (#) temperature coefficient for Na current
    set2_params.Q_10_K = 3; % (#) temperature coefficient for K current 
    
    % Defining the initial Vm, Vmy, n, m, h values (uniformly distributed along the axon)
    set2_params.V_m0 = -58.1124; % (mV) initial condition for membrane potential 
    set2_params.V_my0 = 1.2727; % (mV) initial condition for axon potential in periaxonal space
    set2_params.N_0 = 0.4264; % (#) initial condition for gating variable n
    set2_params.M_0 = 0.1148; % (#) initial condition for gating variable m
    set2_params.H_0 = 0.3548; % (#) initial condition for gating variable h

end
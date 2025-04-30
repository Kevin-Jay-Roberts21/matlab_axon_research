% Set (1) material parameters
% Kevin Roberts
% April 2025

%%%%%%%%%%%%%%%%%%%%%%%%%
% Set (1) Parameter Set %
%%%%%%%%%%%%%%%%%%%%%%%%%
function set1_params = set1_parameters()

    % Defining the material properties
    set1_params.a = 0.55*10^(-4); % (cm) radius in nodal region
    set1_params.a_my = set1_params.a/0.698; % (cm) radius in myelinated region
    set1_params.d_pa = 12.3*10^(-7); % (cm) periaxonal thickness
    set1_params.d_pn = 7.4*10^(-7); % (cm) paranodal thickness
    set1_params.R_i = 0.0712; % (kilo-ohms*cm) intracellular resistivity
    set1_params.R_m = 24.8; % (kilo-ohms*cm^2) specific membrane resistance
    set1_params.C_m = 1.23; % (micro-farads/cm^2) specific membrane capacitance
    set1_params.r_pa = 96.3*10^6; % (kilo-ohms/cm) periaxonal resistivity per unit length
    set1_params.r_pn = 321*10^6; % (kilo-ohms/cm) paranodal resitance per unit length (used in BC since r_bar_pn = r_pn * L_pn) 
    set1_params.R_my = 63.7; % (kilo-ohms*cm^2) specfic myelin resistance
    set1_params.C_my = 0.113; % (micro-fards/cm^2) specific myelin capacitance
    set1_params.R_pa = set1_params.r_pa*pi*set1_params.d_pa*(2*set1_params.a + set1_params.d_pa); % (kilo-ohms*cm) resistivity of the periaxonal space (computed)
    set1_params.R_pn = set1_params.r_pn*pi*set1_params.d_pn*(2*set1_params.a + set1_params.d_pn); % (kilo-ohms*cm) resistivity of the paranodal space (computed)
    set1_params.G_K = 80; % (mS/cm^2) max specific potassium conductance
    set1_params.G_Na = 3000; % (mS/cm^2) max specific sodium conductance 
    set1_params.G_L = 80; % (mS/cm^2) specific leak conductance
    set1_params.E_K = -82; % (mV) Nernst potential for potassium ions
    set1_params.E_Na = 45; % (mV) Nernst potential for sodium ions
    set1_params.E_L = -59.4; % (mV) Nernst potential for leak channels
    set1_params.E_rest = -59.4; % (mV) effective resting nernst potential
    
    % Defining the stimulus value, and time and space interval it is applied to
    set1_params.S_v = 1000; % (mS/cm^2) stimulus value
    set1_params.S_T0 = 1; % (ms) start time of when stimulus is added 
    set1_params.S_T1 = 1.1; % (ms) end time of when stimulus is added 
    set1_params.S_P0 = 0.0001; % (cm) start position of adding the stimulus 
    set1_params.S_P1 = 0.0004; % (cm) end position of adding the stimulus 
    
    % Defining the temperature parameters
    set1_params.T_base = 20; % (C) base temperature
    set1_params.T_actual = 20; % (C) temperature of the squid axon fiber
    set1_params.Q_10_Na = 2.2; % (#) temperature coefficient for Na current
    set1_params.Q_10_K = 3; % (#) temperature coefficient for K current 
    
    % Defining the initial Vm, Vmy, n, m, h values (uniformly distributed along the axon)
    set1_params.V_m0 = -58.1124; % (mV) initial condition for membrane potential 
    set1_params.V_my0 = 1.2727; % (mV) initial condition for axon potential in periaxonal space
    set1_params.N_0 = 0.4264; % (#) initial condition for gating variable n
    set1_params.M_0 = 0.1148; % (#) initial condition for gating variable m
    set1_params.H_0 = 0.3548; % (#) initial condition for gating variable h

end
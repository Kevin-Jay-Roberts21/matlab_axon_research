clear all 
close all
clc

dx = 0.0001; % (cm) space step
dt = 0.01; % (ms) time step 
L_my = 0.0075; % (cm) internodal length
L_n = 0.0005; % (cm) nodal length
L_s = L_n + L_my; % (cm) length of an axon segment
N_s = 7; % (dimless) number of axon segments
L = N_s*L_s; % (cm) total length of axon
T = 5; % (ms) the total time of the experiment
a = 5.5*10^(-5); % (cm) axon radius in nodal region
a_my = 5.623*10^(-5); % (cm) axon radius in myelinated section 
C_m = 0.001; % (ms/(ohms*cm^2)) specific membrane capacitance
C_my = 1.66*10^(-4); % (ms/(ohms*cm^2)) specific myelin capacitance
R_i = 150; % (ohms*cm) intracellular resistivity
R_my = 2.4*10^5; % (ohms*cm^2) specific myelin resistance
R_m = 2.5*10^3; % (ohms*cm^2) specific membrane resistance
G_K = 0.036; % (S/cm^2) max specific potassium conductance
G_Na = 0.12; % (S/cm^2) max specific sodium conductance 
G_L = 0.0003; % (S/cm^2) specific leak conductance
E_K = -82; % (mV) Nernst potential for potassium ions
E_Na = 45; % (mV) Nernst potential for sodium ions
E_L = -59.4; % (mV) Nernst potential for leak channels
V_m0 = ; % (mV) initial condition for membrane potential 
V_my0 = ; % (mV) initial condition for axon potential in periaxonal space
N_0 = ; % (dimless) initial condition for gating variable n
M_0 = ; % (dimless) initial condition for gating variable m
H_0 = ; % (dimless) initial condition for gating variable h

% defining the total number of space steps and time steps
m = L/dx + 1; % m is the number of space steps
n = T/dt; % n is the number of time steps

% defining the alpha_xi and beta_xi functions
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% defining the initial vectors
Vm = V_m0 * ones(1, m);
Vmy = V_my0 * ones(1, m);
N = N_0 * ones(1, m);
M = M_0 * ones(1, m);
H = H_0 * ones(1, m);

% defining some arbitrary constant variables used in discretization
box_1 = a/(2*R_i)*(1 + C_m*a/(C_my*a_my));
box_2 = C_m/(C_my*R_my) - 1/R_m;

for j = 1:n
     
    for i = 2:m-1
        
        
    end
    
    
end
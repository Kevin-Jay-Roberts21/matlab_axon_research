clear all
close all
clc

% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (mS / (ohm*cm^2))
r_l = 30; % membrane resistance (ohms * cm)
a = 0.0025; % axon radius (cm)
d = 5; % axon length (cm)
h = 0.01; % space step (MAY CHANGE LATER)
total_time = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
g_k = 0.036; % (1/(ohm*cm^2))
g_Na = 0.12; % (1/(ohm*cm^2))
g_L = 0.0003; % (1/(ohm*cm^2))
E_k = -77; % (mV)
E_Na = 50; % (mV)
E_L = -54.4; % (mV)

% adding sodium conductance (sitmulus)
S = 0.01;

% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

% can write [U, N, M, H] = axon_simulation_function(...) to get the vectors
% back or just run axon_simulation_function(...) to the just the plots
a = axon_simulation_function(c_m, r_l, a, d, h, total_time, k, g_k, g_Na, g_L, E_k, E_Na, E_L, S, N_0, M_0, H_0, V_initial)

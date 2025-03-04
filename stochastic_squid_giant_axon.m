% Stochastic Hodgkin Huxley model
% This code is meant to add stochasticity to the HH equation through the
% ion channels
% Kevin Roberts
% March 2025

clear all
close all
clc

% Defining all of the material and intrinsic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_m = 1; % membrane capacitance (micro-farads/cm^2)
% R_i = 0.030; % specific intracellular resistivity (kilo-ohms * cm)
% a = 0.025; % axon radius (cm)
G_K = 35; % (mS/cm^2)
G_Na = 40; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % Equilibrium Potential for Potassium Ions (mV)
E_Na = 55; % Equilibrium Potential for Sodium Ions (mV)
E_L = -65; % Equilibrium Potential for Leak Channels (mV)


% Defining the mesh parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 5; % axon length (cm)
% dx = 0.01; % space step (MAY CHANGE LATER)
T = 100; % (ms) Total Time
dt = 0.01; % time step (MAY CHANGE LATER)
% m = L/dx + 1; % total number of space steps
n = T/dt + 1; % total number of time steps

% Defining alpha/beta functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_base = 6.3; % (C) base temperature
T_actual = 6.3; % (C) the temperature of the squid axon
Q_10 = 3; % (dimless) temperature coefficient
phi = Q_10^((T_actual - T_base)/10); % (dimless) temperature scaling factor
alpha_n = @(Vm) phi * 0.02*(Vm - 25)/(1 - exp(-(Vm - 25)/9)); % (1/ms)
alpha_m = @(Vm) phi * 0.182*(Vm + 35)/(1 - exp(-(Vm + 35)/9)); % (1/ms)
alpha_h = @(Vm) phi * 0.25*exp(-(Vm + 90)/12); % (1/ms)
beta_n = @(Vm) phi * -0.002*(Vm - 25)/(1 - exp((Vm - 25)/9)); % (1/ms)
beta_m = @(Vm) phi * -0.124*(Vm + 35)/(1 - exp((Vm + 35)/9)); % (1/ms)
beta_h = @(Vm) phi * 0.25*exp((Vm + 62)/6)/exp((Vm + 90)/12); % (1/ms)

% Stimulus Information
%%%%%%%%%%%%%%%%%%%%%%
S_v = 0; % (in mS/cm^2) % stimulus value
S_T0 = 5; % start time of when stimulus is added (in ms)
S_T1 = 5.1; % end time of when stimulus is added (in ms)
S_P0 = 1; % start position of adding the stimulus (in cm)
S_P1 = 1.1; % end position of adding the stimulus (in cm)

% Initialize Vm, n, m, h
%%%%%%%%%%%%%%%%%%%%%%%%
N_0 = 0.0005632; % probability that potassium gate is open
M_0 = 0.0610; % probability that Sodium activation gate is open
H_0 = 0.5438; % probability that Sodium inactivation gate is open
V_m0 = -60; % (mV) Voltage

% due to the patch clamp, the voltage is the same throughout the axon, i.e.
% at a given time, the voltage is the same at every spatial position.

Vm = zeros(1, n);
N = zeros(1, n);
M = zeros(1, n);
H = zeros(1, n);

Vm(1) = V_m0;
N(1) = N_0;
M(1) = M_0;
H(1) = H_0;

% Starting the time loop
%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(n-1)
    
    xi = randn;

    N(i+1) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
    M(i+1) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
    H(i+1) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
    
    Vm(i+1) = Vm(i) + dt/C_m*sqrt(dt)*xi + dt/C_m*G_Na*(M(i+1))^3*H(i+1)*E_Na + dt/C_m*G_K*(N(i+1))^4*E_K + dt/C_m*G_L*E_L/(1 + dt/C_m*G_Na*(M(i+1))^3*H(i+1) + dt/C_m*G_K*(N(i+1))^4 + dt/C_m*G_L);

end

t = linspace(0, T, n);
plot(t, Vm);




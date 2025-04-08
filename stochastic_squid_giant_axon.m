% Stochastic Hodgkin Huxley model
% This code is meant to add stochasticity to the HH equation through the
% ion channels
% Kevin Roberts
% March 2025

%%
% Reproducing Fig 3 in the resource paper

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
T = 1000; % (ms) Total Time
dt = 0.01; % time step (MAY CHANGE LATER)
% m = L/dx + 1; % total number of space steps
n = T/dt + 1; % total number of time steps
sigma_1 = 30; % the standard deviation of noise (uA/cm^2, microampere per square centimeter) 
sigma_2 = 80;

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
Vm_2 = zeros(1, n);
N_1 = zeros(1, n);
M_1 = zeros(1, n);
H_1 = zeros(1, n);
N_2 = zeros(1, n);
M_2 = zeros(1, n);
H_2 = zeros(1, n);

Vm(1) = V_m0;
Vm_2(1) = V_m0;
N_1(1) = N_0;
M_1(1) = M_0;
H_1(1) = H_0;
N_2(1) = N_0;
M_2(1) = M_0;
H_2(1) = H_0;

% Starting the time loop
%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(n-1)
    
    xi_1 = sigma_1 * sqrt(dt) * randn;
    xi_2 = sigma_2 * sqrt(dt) * randn;

    N_1(i+1) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N_1(i) + dt*alpha_n(Vm(i)));
    M_1(i+1) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M_1(i) + dt*alpha_m(Vm(i)));
    H_1(i+1) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H_1(i) + dt*alpha_h(Vm(i)));
    N_2(i+1) = 1/(1 + dt*alpha_n(Vm_2(i)) + dt*beta_n(Vm_2(i))) * (N_2(i) + dt*alpha_n(Vm_2(i)));
    M_2(i+1) = 1/(1 + dt*alpha_m(Vm_2(i)) + dt*beta_m(Vm_2(i))) * (M_2(i) + dt*alpha_m(Vm_2(i)));
    H_2(i+1) = 1/(1 + dt*alpha_h(Vm_2(i)) + dt*beta_h(Vm_2(i))) * (H_2(i) + dt*alpha_h(Vm_2(i)));

    Vm(i+1) = (Vm(i) + dt/C_m*xi_1 + dt/C_m*G_Na*(M_1(i+1))^3*H_1(i+1)*E_Na + dt/C_m*G_K*(N_1(i+1))^4*E_K + dt/C_m*G_L*E_L)/(1 + dt/C_m*G_Na*(M_1(i+1))^3*H_1(i+1) + dt/C_m*G_K*(N_1(i+1))^4 + dt/C_m*G_L);
    Vm_2(i+1) = (Vm_2(i) + dt/C_m*xi_2 + dt/C_m*G_Na*(M_2(i+1))^3*H_2(i+1)*E_Na + dt/C_m*G_K*(N_2(i+1))^4*E_K + dt/C_m*G_L*E_L)/(1 + dt/C_m*G_Na*(M_2(i+1))^3*H_2(i+1) + dt/C_m*G_K*(N_2(i+1))^4 + dt/C_m*G_L);

end

t = linspace(0, T, n);
plot(t, Vm, 'b-');
hold on;
plot(t, Vm_2, 'r-');

legend('$\sigma = 30 \ \mu $A/cm$^2$', '$\sigma = 80 \ \mu $A/cm$^2$', 'Interpreter', 'latex');

ylabel('Voltage in millivolts.', 'Interpreter', 'latex')
xlabel("Time in milliseconds.")


%%
% Adding N_0, ..., N_4 and M_0, ..., M_3 and H_0, H_1 to the system

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
T = 500; % (ms) Total Time
dt = 0.01; % time step (MAY CHANGE LATER)
% m = L/dx + 1; % total number of space steps
n = T/dt + 1; % total number of time steps
sigma_1 = 10; % the standard deviation of noise (uA/cm^2, microampere per square centimeter) 
N_tot = 11; % number of ion channels open
N_max = 1000; % number of total ion channels

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
N_initial = 0.0005632; % probability that potassium gate is open
M_initial = 0.0610; % probability that Sodium activation gate is open
H_initial = 0.5438; % probability that Sodium inactivation gate is open
V_m0 = -60; % (mV) Voltage

% due to the patch clamp, the voltage is the same throughout the axon, i.e.
% at a given time, the voltage is the same at every spatial position.

Vm = zeros(1, n);
N_0 = zeros(1, n);
N_1 = zeros(1, n);
N_2 = zeros(1, n);
N_3 = zeros(1, n);
N_4 = zeros(1, n);
M_0 = zeros(1, n);
M_1 = zeros(1, n);
M_2 = zeros(1, n);
M_3 = zeros(1, n);
H_0 = zeros(1, n);
H_1 = zeros(1, n);

Vm(1) = V_m0;
N_0(1) = N_initial;
N_1(1) = 0; %N_initial;
N_2(1) = 0; %N_initial;
N_3(1) = 0; %N_initial;
N_4(1) = 0; %N_initial;
M_0(1) = M_initial;
M_1(1) = 0; %M_initial;
M_2(1) = 0; %M_initial;
M_3(1) = 0; %M_initial;
H_0(1) = H_initial;
H_1(1) = 0; % H_initial;

% Starting the time loop
%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(n-1)
    
    xi_1 = sigma_1 * sqrt(dt) * randn;

    N_0(i+1) = (1 + 4*dt*alpha_n(Vm(i)))^(-1)*(N_0(i) + dt*beta_n(Vm(i))*N_1(i));
    N_1(i+1) = (1 + dt*(beta_n(Vm(i)) + 3*alpha_n(Vm(i))))^(-1)*(N_1(i) + dt*(4*alpha_n(Vm(i))*N_0(i) + 2*beta_n(Vm(i))*N_2(i)));
    N_2(i+1) = (1 + dt*(2*beta_n(Vm(i)) + 2*alpha_n(Vm(i))))^(-1)*(N_2(i) + dt*(3*alpha_n(Vm(i))*N_1(i) + 3*beta_n(Vm(i))*N_3(i)));
    N_3(i+1) = (1 + dt*(alpha_n(Vm(i)) + 3*beta_n(Vm(i))))^(-1)*(N_3(i) + dt*(2*alpha_n(Vm(i))*N_2(i) + 4*beta_n(Vm(i))*N_4(i)));
    N_4(i+1) = (1 + 4*dt*beta_n(Vm(i)))^(-1)*(N_4(i) + dt*alpha_n(Vm(i))*N_3(i));
    
    N_sum = N_0(i+1) + N_1(i+1) + N_2(i+1) + N_3(i+1) + N_4(i+1);
    N_0(i+1) = N_0(i+1)/N_sum; 
    N_1(i+1) = N_1(i+1)/N_sum;
    N_2(i+1) = N_2(i+1)/N_sum; 
    N_3(i+1) = N_3(i+1)/N_sum;
    N_4(i+1) = N_4(i+1)/N_sum;

    M_0(i+1) = (1 + 3*dt*alpha_m(Vm(i)))^(-1)*(M_0(i) + dt*beta_m(Vm(i))*M_1(i));
    M_1(i+1) = (1 + dt*(beta_m(Vm(i)) + 2*alpha_m(Vm(i))))^(-1)*(M_1(i) + dt*(3*alpha_m(Vm(i))*M_0(i) + 2*beta_m(Vm(i))*M_2(i)));
    M_2(i+1) = (1 + dt*(alpha_m(Vm(i)) + 2*beta_m(Vm(i))))^(-1)*(M_2(i) + dt*(2*alpha_m(Vm(i))*M_1(i) + 3*beta_m(Vm(i))*M_3(i)));
    M_3(i+1) = (1 + 3*dt*beta_m(Vm(i)))^(-1)*(M_3(i) + dt*alpha_m(Vm(i))*M_2(i));
    
    M_sum = M_0(i+1) + M_1(i+1) + M_2(i+1) + M_3(i+1);
    M_0(i+1) = M_0(i+1)/M_sum; 
    M_1(i+1) = M_1(i+1)/M_sum;
    M_2(i+1) = M_2(i+1)/M_sum; 
    M_3(i+1) = M_3(i+1)/M_sum;

    H_0(i+1) = (1 + dt*alpha_h(Vm(i)))^(-1)*(H_0(i) + dt*(1 - alpha_h(Vm(i)))*H_1(i));
    H_1(i+1) = (1 + dt*(1 - alpha_h(Vm(i))))^(-1)*(H_1(i) + dt*alpha_h(Vm(i))*H_0(i));
    
    H_sum = H_0(i+1) + H_1(i+1);
    H_0(i+1) = H_0(i+1)/H_sum; 
    H_1(i+1) = H_1(i+1)/H_sum;

    G_Na_current = (N_tot / N_max) * G_Na * (M_3(i+1))^3 * H_1(i+1);
    G_K_current = (N_tot / N_max) * G_K * N_4(i+1);
    G_total = G_Na_current + G_K_current + G_L;
    
    I_Na = G_Na_current * E_Na;
    I_K  = G_K_current * E_K;
    I_L  = G_L * E_L;
    
    Vm(i+1) = (Vm(i) + dt/C_m * (xi_1 + I_Na + I_K + I_L)) / (1 + dt/C_m * G_total);

end

t = linspace(0, T, n);
plot(t, Vm, 'b-');
hold on;

legend('$N_{tot} = 11$', 'Interpreter', 'latex');

ylabel('Voltage in millivolts.', 'Interpreter', 'latex')
xlabel("Time in milliseconds.")


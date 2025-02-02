% Computing HH Vm, n, m, and h equilibrium values using Newton-Raphson
% Method
% Kevin Roberts
% Janurary 2025

clear all
close all
clc

syms f1(Vm)

G_K = 36; % (mS/cm^2)
G_Na = 120; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % Equilibrium Potential for Potassium Ions (mV)
E_Na = 50; % Equilibrium Potential for Sodium Ions (mV)
E_L = -54.4; % Equilibrium Potential for Leak Channels (mV)

% SOLVING USING NEWTONS METHOD



% SOLVING USING MATLAB's FSOLVE FUNCTION

% Define the function f(Vm)
f = @(Vm) G_K * (0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10))/(0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80)))^4 * (Vm - E_K) ...
    + G_Na * (0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18)))^3 * ...
    (0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10)))) * (Vm - E_Na) + G_L * (Vm - E_L);

% Set initial guess for Vm
V_m0 = -60; % (mV)

% Use fsolve to find the root of the function (solve for Vm)
options = optimset('Display', 'iter'); % Display iteration details
[Vm_solution, fval] = fsolve(f, V_m0, options);

% Display the solution
disp('Final Vm solution:');
disp(Vm_solution);
disp('Function value at solution (should be close to 0):');
disp(fval);


% Now plugging in the numerically approximated Vm into the probability
% equations:


n_infty = @(Vm) 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10))/(0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80));
m_infty = @(Vm) 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18));
h_infty = @(Vm) 0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10)));

n_infty(Vm_solution)
m_infty(Vm_solution)
h_infty(Vm_solution)


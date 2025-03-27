% Computing HH Vm, n, m, and h equilibrium values using Newton Method and
% fsolve
% Kevin Roberts
% Janurary 2025

clear all
close all
clc

syms f1(Vm)

% % Used for HH Model
G_K = 36; % (mS/cm^2)
G_Na = 500; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % Equilibrium Potential for Potassium Ions (mV)
E_Na = 50; % Equilibrium Potential for Sodium Ions (mV)
E_L = -54.4; % Equilibrium Potential for Leak Channels (mV)

% Used for SC Model
% G_K = 35; % (mS/cm^2)
% G_Na = 40; % (mS/cm^2)
% G_L = 0.3; % (mS/cm^2)
% E_K = -77; % Equilibrium Potential for Potassium Ions (mV)
% E_Na = 55; % Equilibrium Potential for Sodium Ions (mV)
% E_L = -65; % Equilibrium Potential for Leak Channels (mV)

alpha_n = @(Vm) 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)); % (1/ms)
beta_n = @(Vm) 0.125*exp(-(Vm + 65)/80); % (1/ms)
alpha_m = @(Vm) 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)); % (1/ms)
beta_m = @(Vm) 4*exp(-(Vm + 65)/18); % (1/ms)
alpha_h = @(Vm) 0.07*exp(-(Vm + 65)/20); % (1/ms)
beta_h = @(Vm) 1/(1 + exp(-(Vm + 35)/10)); % (1/ms)

% alpha_n = @(Vm) 0.02*(Vm - 25)/(1 - exp(-(Vm - 25)/9)); % (1/ms)
% alpha_m = @(Vm) 0.182*(Vm + 35)/(1 - exp(-(Vm + 35)/9)); % (1/ms)
% alpha_h = @(Vm) 0.25*exp(-(Vm + 90)/12); % (1/ms)
% beta_n = @(Vm) -0.002*(Vm - 25)/(1 - exp((Vm - 25)/9)); % (1/ms)
% beta_m = @(Vm) -0.124*(Vm + 35)/(1 - exp((Vm + 35)/9)); % (1/ms)
% beta_h = @(Vm) 0.25*exp((Vm + 62)/6)/exp((Vm + 90)/12); % (1/ms)


%%
% SOLVING USING MATLAB's FSOLVE FUNCTION

% defining the probability gate function equilibriums
n_infty = @(Vm) alpha_n(Vm)/(alpha_n(Vm) + beta_n(Vm));
m_infty = @(Vm) alpha_m(Vm)/(alpha_m(Vm) + beta_m(Vm));
h_infty = @(Vm) alpha_h(Vm)/(alpha_h(Vm) + beta_h(Vm));


% Define the function f(Vm)
f = @(Vm) G_K*n_infty(Vm)^4*(Vm - E_K) + G_Na*m_infty(Vm)^3*h_infty(Vm)*(Vm - E_Na) + G_L*(Vm - E_L);

% Set initial guess for Vm
V_m0 = -30; % (mV)

% Use fsolve to find the root of the function (solve for Vm)
options = optimset('Display', 'iter'); % Display iteration details
[Vm_solution, fval] = fsolve(f, V_m0, options);

% Display the solution
% disp('Final Vm solution:');
% disp(Vm_solution);
% disp('Function value at solution (should be close to 0):');
% disp(fval);
Vm_solution

% Now plugging in the numerically approximated Vm into the probability
% equations:

n_inf = n_infty(Vm_solution)
m_inf = m_infty(Vm_solution)
h_inf = h_infty(Vm_solution)




function f_of_Vm = f2(Vm, G_K, G_Na, G_L, E_K, E_Na, E_L)

    f_of_Vm = G_K * (0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10))/(0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80)))^4 * (Vm - E_K) ...
    + G_Na * (0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18)))^3 * ...
    (0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10)))) * (Vm - E_Na) + G_L * (Vm - E_L);
end

function derivative_of_f_of_Vm = df2(Vm, G_K, G_Na, G_L, E_K, E_Na, E_L)

    % computed using matlab (this is the result of the computed df1 function)
    derivative_of_f_of_Vm = (36*(Vm/100 + 11/20)^4)/((exp(- Vm/10 - 11/2) - 1)^4*(exp(- Vm/80 - 13/16)/8 - (Vm/100 + 11/20)/(exp(- Vm/10 - 11/2) - 1))^4) + (36*(Vm/100 + 11/20)^3*(Vm + 77))/(25*(exp(- Vm/10 - 11/2) - 1)^4*(exp(- Vm/80 - 13/16)/8 - (Vm/100 + 11/20)/(exp(- Vm/10 - 11/2) - 1))^4) + (144*(Vm/100 + 11/20)^4*(Vm + 77)*(exp(- Vm/80 - 13/16)/640 + 1/(100*(exp(- Vm/10 - 11/2) - 1)) + (exp(- Vm/10 - 11/2)*(Vm/100 + 11/20))/(10*(exp(- Vm/10 - 11/2) - 1)^2)))/((exp(- Vm/10 - 11/2) - 1)^4*(exp(- Vm/80 - 13/16)/8 - (Vm/100 + 11/20)/(exp(- Vm/10 - 11/2) - 1))^5) + (72*exp(- Vm/10 - 11/2)*(Vm/100 + 11/20)^4*(Vm + 77))/(5*(exp(- Vm/10 - 11/2) - 1)^5*(exp(- Vm/80 - 13/16)/8 - (Vm/100 + 11/20)/(exp(- Vm/10 - 11/2) - 1))^4) - (42*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3)/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (63*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^2*(Vm - 50))/(25*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) + (21*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50))/(50*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (126*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50)*((2*exp(- Vm/18 - 65/18))/9 + 1/(10*(exp(- Vm/10 - 4) - 1)) + (exp(- Vm/10 - 4)*(Vm/10 + 4))/(10*(exp(- Vm/10 - 4) - 1)^2)))/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^4) - (42*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*((7*exp(- Vm/20 - 13/4))/2000 - exp(- Vm/10 - 7/2)/(10*(exp(- Vm/10 - 7/2) + 1)^2))*(Vm - 50))/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))^2*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (63*exp(- Vm/10 - 4)*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50))/(25*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^4*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) + 3/10;
end
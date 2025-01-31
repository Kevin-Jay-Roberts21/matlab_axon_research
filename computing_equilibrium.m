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

% 0 = G_K*n_s^4(Vm - E_K) + G_Na*m_s^3*h_s(Vm - E_Na) + G_L(Vm - E_L)

% Thus, f(Vm) = G_K*n_s^4(Vm - E_K) + G_Na*m_s^3*h_s(Vm - E_Na) + G_L(Vm - E_L)

% where n_s = alpha_n(Vm)/(alpha_n(Vm) + beta_n(Vm))  which implies...
%       n_s = 0.01*(Vm + 55)/((1 - exp(-(Vm + 55)/10))/0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80))
%       
%       m_s = alpha_m(Vm)/(alpha_m(Vm) + beta_m(Vm))  which implies...
%       m_s = 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18))
%       
%       h_s = alpha_h(Vm)/(alpha_h(Vm) + beta_h(Vm))  which implies...
%       h_s = 0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10)))


f1(Vm) = G_K*(0.01*(Vm + 55)/((1 - exp(-(Vm + 55)/10))/0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80)))^4*(Vm - E_K) + G_Na*(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18)))^3*(0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10))))*(Vm - E_Na) + G_L*(Vm - E_L);
df1(Vm) = diff(f1, Vm, 1); % the first derivative of f1(Vm) w.r.t. Vm

function f_of_Vm = f2(Vm, G_K, G_Na, G_L, E_K, E_Na, E_L)

    f_of_Vm = G_K*(0.01*(Vm + 55)/((1 - exp(-(Vm + 55)/10))/0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)) + 0.125*exp(-(Vm + 65)/80)))^4*(Vm - E_K) + G_Na*(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10))/(0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)) + 4*exp(-(Vm + 65)/18)))^3*(0.07*exp(-(Vm + 65)/20)/(0.07*exp(-(Vm + 65)/20) + 1/(1 + exp(-(Vm + 35)/10))))*(Vm - E_Na) + G_L*(Vm - E_L);
end

function derivative_of_f_of_Vm = df2(Vm, G_K, G_Na, G_L, E_K, E_Na, E_L)

    % computed using matlab (this is the result of the computed df1 function)
    derivative_of_f_of_Vm = (36*(Vm/100 + 11/20)^4)/(exp(- Vm/80 - 13/16)/8 + ((100*exp(- Vm/10 - 11/2) - 100)*(Vm + 55))/(exp(- Vm/10 - 11/2) - 1))^4 + (36*(Vm/100 + 11/20)^3*(Vm + 77))/(25*(exp(- Vm/80 - 13/16)/8 + ((100*exp(- Vm/10 - 11/2) - 100)*(Vm + 55))/(exp(- Vm/10 - 11/2) - 1))^4) + (144*(Vm/100 + 11/20)^4*(Vm + 77)*(exp(- Vm/80 - 13/16)/640 - (100*exp(- Vm/10 - 11/2) - 100)/(exp(- Vm/10 - 11/2) - 1) + (10*exp(- Vm/10 - 11/2)*(Vm + 55))/(exp(- Vm/10 - 11/2) - 1) - (exp(- Vm/10 - 11/2)*(100*exp(- Vm/10 - 11/2) - 100)*(Vm + 55))/(10*(exp(- Vm/10 - 11/2) - 1)^2)))/(exp(- Vm/80 - 13/16)/8 + ((100*exp(- Vm/10 - 11/2) - 100)*(Vm + 55))/(exp(- Vm/10 - 11/2) - 1))^5 - (42*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3)/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (63*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^2*(Vm - 50))/(25*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) + (21*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50))/(50*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (126*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50)*((2*exp(- Vm/18 - 65/18))/9 + 1/(10*(exp(- Vm/10 - 4) - 1)) + (exp(- Vm/10 - 4)*(Vm/10 + 4))/(10*(exp(- Vm/10 - 4) - 1)^2)))/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^4) - (42*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*((7*exp(- Vm/20 - 13/4))/2000 - exp(- Vm/10 - 7/2)/(10*(exp(- Vm/10 - 7/2) + 1)^2))*(Vm - 50))/(5*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))^2*(exp(- Vm/10 - 4) - 1)^3*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) - (63*exp(- Vm/10 - 4)*exp(- Vm/20 - 13/4)*(Vm/10 + 4)^3*(Vm - 50))/(25*((7*exp(- Vm/20 - 13/4))/100 + 1/(exp(- Vm/10 - 7/2) + 1))*(exp(- Vm/10 - 4) - 1)^4*(4*exp(- Vm/18 - 65/18) - (Vm/10 + 4)/(exp(- Vm/10 - 4) - 1))^3) + 3/10;
end

% initial V_m guess
V_m0 = -50; % (mV)
V_m = V_m0;

% Preallocate array to store Vm values
Vm_values = zeros(1, 100);

% damping factor
alpha = 0.01;

for i = 1:10000

    newVm = V_m - alpha * f2(V_m, G_K, G_Na, G_L, E_K, E_Na, E_L)/df2(V_m, G_K, G_Na, G_L, E_K, E_Na, E_L);
    
    % Store the new Vm value
    Vm_values(i) = newVm;

    % update
    V_m = newVm;

     % Plot every 100th Vm value to avoid overloading the plot
    plot(1:i, Vm_values(1:i));
    xlabel('Iteration');
    ylabel('V_m (mV)');
        
end

final_Vm = newVm


%% 
% OR using matlab's built in Newtons method:

% Define the function f(Vm)
f = @(Vm) G_K * (0.01 * (Vm + 55) / ((1 - exp(-(Vm + 55)/10)) / (0.01 * (Vm + 55)) / (1 - exp(-(Vm + 55)/10)) + 0.125 * exp(-(Vm + 65)/80)))^4 * (Vm - E_K) ...
    + G_Na * (0.1 * (Vm + 40) / (1 - exp(-(Vm + 40)/10)) / (0.1 * (Vm + 40) / (1 - exp(-(Vm + 40)/10)) + 4 * exp(-(Vm + 65)/18)))^3 * ...
    (0.07 * exp(-(Vm + 65)/20) / (0.07 * exp(-(Vm + 65)/20) + 1 / (1 + exp(-(Vm + 35)/10)))) * (Vm - E_Na) + G_L * (Vm - E_L);

% Set initial guess for Vm
V_m0 = -50; % (mV)

% Use fsolve to find the root of the function (solve for Vm)
options = optimset('Display', 'iter'); % Display iteration details
[Vm_solution, fval] = fsolve(f, V_m0, options);

% Display the solution
disp('Final Vm solution:');
disp(Vm_solution);
disp('Function value at solution (should be close to 0):');
disp(fval);

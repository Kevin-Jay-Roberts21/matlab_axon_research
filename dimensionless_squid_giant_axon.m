clear all
close all
clc

% this is the nondimensionalized discretization of the HH PDE for squid
% giant axons


% Defining all of the material and intrinsic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_m = 1; % membrane capacitance (micro-farads/cm^2)
R_i = 0.03; % specific intracellular resistivity (kilo-ohms * cm)
a = 0.025; % axon radius (cm)
G_K = 36; % (mS/cm^2)
G_Na = 120; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % Equilibrium Potential for Potassium Ions (mV)
E_Na = 50; % Equilibrium Potential for Sodium Ions (mV)
E_L = -54.4; % Equilibrium Potential for Leak Channels (mV)

% Defining the mesh parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 5; % axon length (cm)
dx = 0.01; % space step (MAY CHANGE LATER)
T = 20; % (ms) Total Time
dt = 0.01; % time step (MAY CHANGE LATER)

% Defining alpha/beta functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_base = 6.3; % (C) base temperature
T_actual = 6.3; % (C) the temperature of the squid axon
Q_10 = 3; % (dimless) temperature coefficient
phi = Q_10^((T_actual - T_base)/10); % (dimless) temperature scaling factor
alpha_n = @(Vm) phi * 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)); % (1/ms)
beta_n = @(Vm) phi * 0.125*exp(-(Vm + 65)/80); % (1/ms)
alpha_m = @(Vm) phi * 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)); % (1/ms)
beta_m = @(Vm) phi * 4*exp(-(Vm + 65)/18); % (1/ms)
alpha_h = @(Vm) phi * 0.07*exp(-(Vm + 65)/20); % (1/ms)
beta_h = @(Vm) phi * 1/(1 + exp(-(Vm + 35)/10)); % (1/ms)

% Stimulus Information
%%%%%%%%%%%%%%%%%%%%%%
S_v = 10; % (in mS/cm^2) % stimulus value
S_T0 = 5; % start time of when stimulus is added (in ms)
S_T1 = 5.1; % end time of when stimulus is added (in ms)
S_P0 = 1; % start position of adding the stimulus (in cm)
S_P1 = 1.1; % end position of adding the stimulus (in cm)

% Defining the Dimensionless Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% characteristic scaling constants
t_c = 1; % (in ms)
x_c = 5; % (in cm)
V_c = 25; % (k_B*T/e = 25mV)

% material params converted to dimensionless
G_K_tilde = t_c*G_K/C_m;
G_Na_tilde = t_c*G_Na/C_m;
G_L_tilde = t_c*G_L/C_m;
E_K_tilde = E_K/V_c;
E_Na_tilde = E_Na/V_c;
E_L_tilde = E_L/V_c;

% mesh params converted to dimensionless
T_tilde = T/t_c;
L_tilde = L/x_c;

% alpha/beta converted to dimensionless
alpha_n_tilde = @(V) t_c*alpha_n(V*V_c);
beta_n_tilde = @(V) t_c*beta_n(V*V_c);
alpha_m_tilde = @(V) t_c*alpha_m(V*V_c);
beta_m_tilde = @(V) t_c*beta_m(V*V_c);
alpha_h_tilde = @(V) t_c*alpha_h(V*V_c);
beta_h_tilde = @(V) t_c*beta_h(V*V_c);

% stimulus info converted to dimensionless
S_T0_tilde = S_T0/t_c;
S_T1_tilde = S_T1/t_c;
S_P0_tilde = S_P0/x_c; 
S_P1_tilde = S_P1/x_c;
S_v_tilde = t_c*S_v/C_m;

% m and n must be defined using dimless L and T
m = L_tilde/dx + 1; % total number of space steps
n = T_tilde/dt + 1; % total number of time steps

% Initialize Vm, n, m, h
%%%%%%%%%%%%%%%%%%%%%%%%
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_m0 = -64.9997/V_c; % dimenionless voltage
Vm = V_m0 * ones(1, m);
N = N_0 * ones(1, m);
M = M_0 * ones(1, m);
H = H_0 * ones(1, m);

% defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

% Initialize matrix A and vectors b and f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(m);
b = zeros(m, 1);
f = zeros(m, 1);

% Starting the time loop
%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)

    % setting newN, newM, newH vectors (this is a new i, different from the above for loop)
    for i = 1:m
        newN(i) = 1/(1 + dt*alpha_n_tilde(Vm(i)) + dt*beta_n_tilde(Vm(i))) * (N(i) + dt*alpha_n_tilde(Vm(i)));
        newM(i) = 1/(1 + dt*alpha_m_tilde(Vm(i)) + dt*beta_m_tilde(Vm(i))) * (M(i) + dt*alpha_m_tilde(Vm(i)));
        newH(i) = 1/(1 + dt*alpha_h_tilde(Vm(i)) + dt*beta_h_tilde(Vm(i))) * (H(i) + dt*alpha_h_tilde(Vm(i)));
    end

    % Edit the next Vm, N, M and H (redefining Vm, N, M and H vectors)
    N = newN;
    M = newM;
    H = newH;
    
    % i is the space step
    for i = 1:m
        
        a1 = -dt*a/(2*R_i*C_m*dx^2) * t_c/x_c^2;
        a2 = 1 + dt*a/(R_i*C_m*dx^2) * t_c/x_c^2 + dt*G_K_tilde*N(i)^4 + dt*G_Na_tilde*M(i)^3*H(i) + dt*G_L_tilde;
        a3 = -dt*a/(2*R_i*C_m*dx^2) * t_c/x_c^2;
        a4 = 1; 
        a5 = dt*G_K_tilde*N(i)^4*E_K_tilde + dt*G_Na_tilde*M(i)^3*H(i)*E_Na_tilde + dt*G_L_tilde*E_L_tilde;
        
        % adding stimulus temporally and spatially:
        if (j*dt >= S_T0_tilde && j*dt <= S_T1_tilde) && (i*dx >= S_P0_tilde && i*dx <= S_P1_tilde)
            a2 = 1 + dt*a/(R_i*C_m*dx^2) * t_c/x_c^2 + dt*G_K_tilde*N(i)^4 + dt*(G_Na_tilde*M(i)^3*H(i)+ S_v_tilde) + dt*G_L_tilde;
            a5 = dt*G_K_tilde*N(i)^4*E_K_tilde + dt*(G_Na_tilde*M(i)^3*H(i) + S_v_tilde)*E_Na_tilde + dt*G_L_tilde*E_L_tilde;
        end

        % constructing the A matrix and b and f vectors
        if i == 1
            A(1, 1) = 1;
            A(1, 2) = -1;
            b(1) = 0;
            f(1) = 0;
        elseif i == m
            A(m, m-1) = -1;
            A(m, m) = 1;
            b(m) = 0;
            f(m) = 0;
        else
            A(i, i-1) = a3;
            A(i, i) = a2;
            A(i, i+1) = a1;
            b(i) = a4*Vm(i); 
            f(i) = a5;
        end
    end

    % setting newVm
    newVm = transpose(A\(b+f));
    Vm = newVm;

    % adding the newly defined vectors to the 'all' matrices
    Vm_all(j+1,:) = Vm;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

end

% now pick a position to plot all of the voltages
position1 = 0.5/x_c; % (dimless)
position2 = 1/x_c; % (dimless)
position3 = 1.5/x_c; % (dimless)
position4 = 2/x_c; % (dimless)
position5 = 2.5/x_c; % (dimless)
position6 = 3/x_c; % (dimless)
position7 = 3.5/x_c; % (dimless)

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7];

% Times to observe the voltage along the axon
time1 = 5/t_c; % in ms
time2 = 6/t_c; % in ms
time3 = 8/t_c; % in ms
time4 = 10/t_c; % in ms
time5 = 11/t_c; % in ms
time6 = 12/t_c; % in ms
time7 = 13/t_c; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];

% plotting Voltage vs Axon length             
figure(1)
t1 = linspace(0, L_tilde, m);
plot(t1, Vm_all(round(time1/dt),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Vm_all(round(list_of_times(i)/dt),:))
end

% describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('$\\tilde{V_m}$ at $\\tilde{t} = %g$', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex')
ylabel("Dimensionless Voltage $\tilde{V_m}$.", 'Interpreter','latex')
xlabel("Dimensionless length $\tilde{d}$ of the axon.", 'Interpreter','latex')

% plotting Voltage vs Time
figure(2)
t2 = linspace(0, T_tilde, n); % FULL MATRIX
% t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Vm_all(:,round(position1/dx)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, Vm_all(:,round(list_of_positions(i)/dx)))
end

% describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('$\\tilde{V_m}$ at $\\tilde{x} = %g$', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel("Dimensionless Voltage $\tilde{V_m}$.", 'Interpreter','latex')
xlabel("Dimensionless Time $\tilde{T}$.", 'Interpreter','latex')

%plotting N, M, H probability vs time (at a certain position)
figure(3)
plot(t2, N_all(:,round(position3/dx)))
hold on
plot(t2, M_all(:,round(position3/dx)))
hold on
plot(t2, H_all(:,round(position3/dx)))
legendStrings3 = {
    sprintf('N at $\\tilde{x} = %g$', position3), ...
    sprintf('M at $\\tilde{x} = %g$', position3), ...
    sprintf('H at $\\tilde{x} = %g$', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Dimensionless Time $\tilde{T}$.", 'Interpreter','latex')

% save('dimless_stim_0.003912.mat');
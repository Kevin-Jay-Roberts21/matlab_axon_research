clear all
close all
clc

% this is the nondimensionalized discretization of the HH PDE for squid
% giant axons


% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (ohms * cm)
a = 0.0025; % axon radius (cm)
L = 5; % axon length (cm)
T = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
g_k = 0.036; % (1/(ohm*cm^2))
g_Na = 0.12; % (1/(ohm*cm^2))
g_L = 0.0003; % (1/(ohm*cm^2))
E_k = -77; % (mV)
E_Na = 50; % (mV)
E_L = -54.4; % (mV)

% The following are functions of voltage
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% defining the characteristic scaling constants
t_c = 1; % (in ms)
x_c = 5; % (in cm)
V_c = 25; % (k_B*T/e = 25mV)

% Defining all of the dimensionless parameters
gamma = t_c*a/(2*r_l*c_m*x_c^2);
g_k_tilde = t_c*g_k/c_m;
g_Na_tilde = t_c*g_Na/c_m;
g_L_tilde = t_c*g_L/c_m;
E_k_tilde = E_k/V_c;
E_Na_tilde = E_Na/V_c;
E_L_tilde = E_L/V_c;
alpha_n_tilde = @(V) t_c*alpha_n(V*V_c);
beta_n_tilde = @(V) t_c*beta_n(V*V_c);
alpha_m_tilde = @(V) t_c*alpha_m(V*V_c);
beta_m_tilde = @(V) t_c*beta_m(V*V_c);
alpha_h_tilde = @(V) t_c*alpha_h(V*V_c);
beta_h_tilde = @(V) t_c*beta_h(V*V_c);
T_tilde = T/t_c; % since t_c = 1, nothing changes, but this is still good to do if t_c is not 1.
L_tilde = L/x_c;

h = 0.002; % space step
k = 0.01; % time step

% adding sodium conductance (stimulus)
% recall that we are in terms of tau and chi now, namely, if I had the time
% interval as T0-T1 before, now it should be: T0_tau = T0/t_c -
% T1_tau = T1/t_c AND for the space interval if it was P0-P1, now it should
% be P0_chi = P0/x_c - P1_chi = P1/x_c.
S = 0.003912; % (in 1/(ohm*cm^2))
T0 = 5; % start time of when stimulus is added (in ms)
T1 = 5.1; % end time of when stimulus is added (in ms)
P0 = 1; % position of adding the stimulus (in cm)
P1 = 1.1;

% converted to dimensionless
T0_tilde = T0/t_c;
T1_tilde = T1/t_c;
P0_tilde = P0/x_c; 
P1_tilde = P1/x_c;
S_tilde = t_c*S/c_m;


% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997/V_c; % (mV) Voltage (eq: -64.9997) new equilibrium is: -2.59999 

% number of columns of the matrices (length of axon divided by space step)
m = L_tilde/h + 1; 

% initial vectors
U = zeros(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);

% A and b matrix and vector (will be solved in Ax = b style)
A = zeros(m);
b = zeros(m, 1);

% setting the initial U H M and N conditions in the vectors:
U(1,:) = V_initial;
N(1,:) = N_0;
M(1,:) = M_0;
H(1,:) = H_0;

% Final matri with voltage of axon at every time and space step
Uall(1,:) = U;

% Final matrix with probabilities at every time and space step
Nall(1,:) = N;
Mall(1,:) = M; 
Hall(1,:) = H;

% defining new time vector (each time vector has m elements)
newU = zeros(1, m);
newN = zeros(1, m);
newM = zeros(1, m);
newH = zeros(1, m);

% number of rows in final matrices
n = T_tilde/k;

% j is the time step
for j = 1:(n-1)

    % setting newN, newM, newH vectors (this is a new i, different from the above for loop)
    for i = 1:m
        newN(i) = 1/(1/k + alpha_n_tilde(U(i)) + beta_n_tilde(U(i))) * (N(i)/k + alpha_n_tilde(U(i)));
        newM(i) = 1/(1/k + alpha_m_tilde(U(i)) + beta_m_tilde(U(i))) * (M(i)/k + alpha_m_tilde(U(i)));
        newH(i) = 1/(1/k + alpha_h_tilde(U(i)) + beta_h_tilde(U(i))) * (H(i)/k + alpha_h_tilde(U(i)));
    end

    % Edit the next U, N, M and H (redefining U, N, M and H vectors)
    N = newN;
    M = newM;
    H = newH;
    
    % i is the space step
    for i = 1:m
        
        % defining coefficients
        a1 = -k*gamma/h^2;
        a2 = 1 + 2*k*gamma/h^2 + k*g_k_tilde*N(i)^4 + k*g_Na_tilde*M(i)^3*H(i) + k*g_L_tilde;
        a3 = -k*gamma/h^2;
        a4 = 1; 
        a5 = k*g_k_tilde*N(i)^4*E_k_tilde + k*g_Na_tilde*M(i)^3*H(i)*E_Na_tilde + k*g_L_tilde*E_L_tilde;

        % % adding the stimulus temporally only: (T0 - T1)
        % if j*k >= T0_tilde && j*k <= T1_tilde
        %     a2 = 1 + 2*k*gamma/h^2 + k*g_k_tilde*N(i)^4 + k*(g_Na_tilde*M(i)^3*H(i)+ S_tilde) + k*g_L_tilde;
        %     a5 = k*g_k_tilde*N(i)^4*E_k_tilde + k*(g_Na_tilde*M(i)^3*H(i) + S_tilde)*E_Na_tilde + k*g_L_tilde*E_L_tilde;
        % end
        
        % adding the stimulus spatially only: (P0 - P1)
        % if i*h >= P0_tilde && i*h <= P1_tilde 
        %     a2 = 1 + 2*k*gamma/h^2 + k*g_k_tilde*N(i)^4 + k*(g_Na_tilde*M(i)^3*H(i)+ S_tilde) + k*g_L_tilde;
        %     a5 = k*g_k_tilde*N(i)^4*E_k_tilde + k*(g_Na_tilde*M(i)^3*H(i) + S_tilde)*E_Na_tilde + k*g_L_tilde*E_L_tilde;
        % end
        
        % adding stimulus temporally and spatially:
        if (j*k >= T0_tilde && j*k <= T1_tilde) && (i*h >= P0_tilde && i*h <= P1_tilde)
            a2 = 1 + 2*k*gamma/h^2 + k*g_k_tilde*N(i)^4 + k*(g_Na_tilde*M(i)^3*H(i)+ S_tilde) + k*g_L_tilde;
            a5 = k*g_k_tilde*N(i)^4*E_k_tilde + k*(g_Na_tilde*M(i)^3*H(i) + S_tilde)*E_Na_tilde + k*g_L_tilde*E_L_tilde;
        end

        % constructing the A and b matrix and vector
        if i == 1
            A(1, 1) = 1/h;
            A(1, 2) = -1/h;
            b(1, 1) = 0;
        elseif i == m
            A(m, m-1) = -1/h;
            A(m, m) = 1/h;
            b(m, 1) = 0;
        else
            A(i, i-1) = a1;
            A(i, i) = a2;
            A(i, i+1) = a3;
            b(i, 1) = a4*U(i) + a5; 
        end
    end

    % setting newU (the solution from Ax = b)
    newU = transpose(A\b);
    U = newU;

    % adding the newly defined vectors to the 'all' matrices
    Uall(j+1,:) = U;
    Nall(j+1,:) = N;
    Mall(j+1,:) = M;
    Hall(j+1,:) = H;

% USE FOR GRABBING AT EVERY 50th ITERATION
%     if mod(j, 50) == 0
%         Uall(round(j/50)+1,:) = U;
%         Nall(round(j/50)+1,:) = N;
%         Mall(round(j/50)+1,:) = M;
%         Hall(round(j/50)+1,:) = H;
%     end
% USE FOR GRABBING AT EVERY 100th ITERATION
%     if mod(j, 100) == 0
%         Uall(round(j/100)+1,:) = U;
%         Nall(round(j/100)+1,:) = N;
%         Mall(round(j/100)+1,:) = M;
%         Hall(round(j/100)+1,:) = H;
%     end

end

% now pick a position to plot all of the voltages
position1 = 0.1; % in cm
position2 = 0.2; % in cm
position3 = 0.3; % in cm
position4 = 0.4; % in cm
position5 = 0.5; % in cm
position6 = 0.6; % in cm
position7 = 0.7; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7];

% Times to observe the voltage along the axon
time1 = 5; % in ms
time2 = 6; % in ms
time3 = 8; % in ms
time4 = 10; % in ms
time5 = 11; % in ms
time6 = 12; % in ms
time7 = 13; % in ms

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
plot(t1, Uall(round(time1/k),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Uall(round(list_of_times(i)/k),:))
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
t2 = linspace(0, T, n); % FULL MATRIX
% t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Uall(:,round(position1/h)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, Uall(:,round(list_of_positions(i)/h)))
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
plot(t2, Nall(:,round(position3/h)))
hold on
plot(t2, Mall(:,round(position3/h)))
hold on
plot(t2, Hall(:,round(position3/h)))
legendStrings3 = {
    sprintf('N at $\\tilde{x} = %g$', position3), ...
    sprintf('M at $\\tilde{x} = %g$', position3), ...
    sprintf('H at $\\tilde{x} = %g$', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Dimensionless Time $\tilde{T}$.", 'Interpreter','latex')

% save('dimless_stim_0.003912.mat');
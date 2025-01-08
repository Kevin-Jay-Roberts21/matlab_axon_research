% "Myelinated" Hodgkin Huxley Equation Using Finite Difference Methods and
% spatial functions for radius, capacitance, Na conductance and K
% conductance
% Kevin Robert
% December 2024

clear all
close all
clc

% Defining all of the material and intrinsic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_m = 1; % membrane capacitance (micro-farads/cm^2)
C_my = 1; % myelin capacitance (micro-farads/cm^2)
R_i = 0.03; % specific intracellular resistivity (or axoplasmic resistivity) (kilo-ohms * cm)
a = 0.00005; % (cm)
a_my = 0.00005; % (cm)
G_K_nodal = 36; % (mS/cm^2)
G_K_internodal = 0; % (mS/cm^2)
G_Na_nodal = 500; % (mS/cm^2)
G_Na_internodal = 0; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % (mV)
E_Na = 50; % (mV)
E_L = -54.4; % (mV)

% Defining the Mesh Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.0005; % space step (cm)
dt = 0.01; % time step (MAY CHANGE LATER)
L_my = 0.0115; % (in cm)
L_n = 0.0005; % (in cm)
n_s = 15; % total number of axon segments
L = (L_my*(n_s-1)) + (L_n*(n_s-1)) + L_n; % axon length (in cm)
T = 7; % (ms) Total Time
m = round(L/dx) + 1; % Total number of space steps
n = T/dt + 1; % Total number of time steps

% Defining the alpha/beta functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
S_v = 10; % (in mS/cm^2)
S_T0 = 0.1; % start time of when stimulus is added (in ms)
S_T1 = 0.2; % end time of when stimulus is added (in ms)
S_P0 = 0.0120; % position of adding the stimulus (in cm)
S_P1 = 0.0125; % ending position of adding the stimulus (in cm or *10^4 in um)

% Creating the nodal and internodal functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating a list of nodal regions [[start_pos1, end_pos1], [start_pos2, end_pos2], ...]
% reasoning for this here is to use nodal_regions to describe the 4 spacial
% functions (radius, capacitance, K conductance and Na conductance)
nodal_regions = [];
for i = 0:(n_s-1)
    nodal_regions(:,i+1) = [(i*L_n)+(i*L_my), (i*L_n)+(i*L_my)+L_n];
end

% The following are functions of space starting each of the functions with 
% a condition. This condition being the values of each variable if the position 
% is in the first nodal region
a_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
C_m_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
G_Na_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
G_K_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));

% Creating a for loop that will turn the a, c_m, g_k and g_Na into piecewise
% functions, having different values in nodal regions vs internodal regions 
for i = 2:n_s
    a_fctn = @(x) a_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    C_m_fctn = @(x) C_m_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    G_Na_fctn = @(x) G_Na_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    G_K_fctn = @(x) G_K_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
end

% Finally, we make a, g_Na, and g_k certain values if the conditions above
% are met or not the first element of the summation is if x is in a nodal
% region, the second element is if the x is in a myelinated region.
a_fctn = @(x) a*a_fctn(x) + a_my*(a_fctn(x)==0); 
C_m_fctn = @(x) C_m*C_m_fctn(x) + C_my*(C_m_fctn(x)==0);
G_K_fctn = @(x) G_K_nodal*G_K_fctn(x) + G_K_internodal*(G_K_fctn(x)==0); 
G_Na_fctn = @(x) G_Na_nodal*G_Na_fctn(x) + G_Na_internodal*(G_Na_fctn(x)==0); 

% Initialize Vm, n, m, and h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
% M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
% H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
% V_m0 = -90; % (mV) Voltage (eq: -64.9997)
N_0 = 0.5065; % probability that potassium gate is open 
M_0 = 0.1919; % probability that Sodium activation gate is open 
H_0 = 0.2128; % probability that Sodium inactivation gate is open
V_m0 = -52.9326; % (mV) Voltage 
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
        newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
        newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
        newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
    end

    % Edit the next Vm, N, M and H (redefining Vm, N, M and H vectors)
    N = newN;
    M = newM;
    H = newH;
    
    % i is the space step
    for i = 1:m

        a1 = -dt*a_fctn(i*dx)/(2*R_i*C_m_fctn(i*dx)*dx^2);
        a2 = 1 + dt*a_fctn(i*dx)/(R_i*C_m_fctn(i*dx)*dx^2) + dt*G_K_fctn(i*dx)*(N(i)^4)/C_m_fctn(i*dx) + dt*G_Na_fctn(i*dx)*(M(i)^3)*H(i)/C_m_fctn(i*dx) + dt*G_L/C_m_fctn(i*dx);
        a3 = -dt*a_fctn(i*dx)/(2*R_i*C_m_fctn(i*dx)*dx^2); 
        a4 = 1; 
        a5 = dt*G_K_fctn(i*dx)*(N(i)^4)*E_K/C_m_fctn(i*dx) + dt*G_Na_fctn(i*dx)*(M(i)^3)*H(i)*E_Na/C_m_fctn(i*dx) + dt*G_L*E_L/C_m_fctn(i*dx);
        
        % adding stimulus temporally and spatially:
        if (j*dt >= S_T0 && j*dt <= S_T1) && (i*dx >= S_P0 && i*dx <= S_P1)
            a2 = 1 + dt*a_fctn(i*dx)/(R_i*C_m_fctn(i*dx)*dx^2) + dt*G_K_fctn(i*dx)*(N(i)^4)/C_m_fctn(i*dx) + dt*(G_Na_fctn(i*dx)*(M(i)^3)*H(i) + S_v)/C_m_fctn(i*dx) + dt*G_L/C_m_fctn(i*dx);
            a5 = dt*G_K_fctn(i*dx)*(N(i)^4)*E_K/C_m_fctn(i*dx) + dt*(G_Na_fctn(i*dx)*(M(i)^3)*H(i) + S_v)*E_Na/C_m_fctn(i*dx) + dt*G_L*E_L/C_m_fctn(i*dx);
        end

        % constructing the A and b matrix and vector
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

    % setting newU (the solution from Ax = b))
    newVm = transpose(A\(b+f));
    Vm = newVm;

    % adding the newly defined vectors to the 'all' matrices
    Vm_all(j+1,:) = Vm;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

end


% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = 0.01; % in cm 
position2 = 0.02;
position3 = 0.03; 
position4 = 0.04; 
position5 = 0.05; 
position6 = 0.06; 
position7 = 0.09; 
position8 = 0.12;
position9 = 0.14;
position10 = 0.16;

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7
                     position8
                     position9
                     position10];

% Times to observe the voltage along the axon
time1 = 2; % in ms
time2 = 2.1; % in ms
time3 = 3; % in ms
time4 = 3; % in ms
time5 = 4; % in ms
time6 = 5; % in ms
time7 = 6; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];

% plotting Voltage vs Axon length
figure(1)
t1 = linspace(0, L, m);
plot(t1, Vm_all(round(time1/dt),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Vm_all(round(list_of_times(i)/dt),:))
end

% describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('$V_m$ at t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex')
ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
xlabel("Length of the axon in cm.")

% plotting Voltage vs Time
figure(2)
t2 = linspace(0, T, n); % FULL MATRIX
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
    legendStrings2{end+1} = sprintf('$V_m$ at x = %g cm', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
xlabel("Time in milliseconds.")

% plotting N, M, H probability vs time (at a certain position)
figure(3)
plot(t2, N_all(:,round(position1/dx)))
hold on
plot(t2, M_all(:,round(position1/dx)))
hold on
plot(t2, H_all(:,round(position1/dx)))
legendStrings3 = {
    sprintf('n at x = %g cm', position1), ...
    sprintf('m at x = %g cm', position1), ...
    sprintf('h at x = %g cm', position1)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")

% save('myelin_stim_0.5_radius_0.00005_starting_at_equil_15_node_h_.00001.mat');
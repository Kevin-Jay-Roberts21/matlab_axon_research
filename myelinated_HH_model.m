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
a = 0.025; % (cm)
a_my = 0.025; % (cm)
G_K_nodal = 36; % (mS/cm^2)
G_K_internodal = 36; % (mS/cm^2)
G_Na_nodal = 500; % (mS/cm^2)
G_Na_internodal = 500; % (mS/cm^2)
G_L = 0.3; % (mS/cm^2)
E_K = -77; % (mV)
E_Na = 50; % (mV)
E_L = -54.4; % (mV)

% Defining the Mesh Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.00005; % space step (cm)
dt = 0.01; % time step (MAY CHANGE LATER)
L_my = 0.0075; % (in cm)
L_n = 0.0005; % (in cm)
L_s = L_n + L_my; % (cm) length of an axon segment
n_s = 50; % total number of axon segments
L = n_s*L_s; % axon length (in cm)
T = 10; % (ms) Total Time
N_n = round(L_n/dx); % number of space steps in a nodal region
N_my = round(L_my/dx); % number of space steps in an internodal region
N_s = N_n + N_my; % number of space steps in an entire axon segement
m = N_s*n_s + 1; % total number of space steps
n = T/dt + 1; % n is the number of time steps

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
S_v = 100; % (in mS/cm^2)
S_T0 = 0.1; % start time of when stimulus is added (in ms)
S_T1 = 0.2; % end time of when stimulus is added (in ms)
S_P0 = 0.0001; % position of adding the stimulus (in cm)
S_P1 = 0.0004; % ending position of adding the stimulus (in cm or *10^4 in um)

% Creating the nodal and internodal functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating a list of nodal regions [[start_pos1, end_pos1], [start_pos2, end_pos2], ...]
% reasoning for this here is to use nodal_regions to describe the 4 spacial
% functions (radius, capacitance, K conductance and Na conductance)
% nodal_regions = [];
% for i = 0:(n_s-1)
%     nodal_regions(:,i+1) = [(i*L_n)+(i*L_my), (i*L_n)+(i*L_my)+L_n];
% end

% The following are functions of space starting each of the functions with 
% a condition. This condition being the values of each variable if the position 
% is in the first nodal region
% a_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
% C_m_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
% G_Na_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
% G_K_fctn = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));

% Creating a for loop that will turn the a, c_m, g_k and g_Na into piecewise
% functions, having different values in nodal regions vs internodal regions 
% for i = 2:n_s
%     a_fctn = @(x) a_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
%     C_m_fctn = @(x) C_m_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
%     G_Na_fctn = @(x) G_Na_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
%     G_K_fctn = @(x) G_K_fctn(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
% end

% Finally, we make a, g_Na, and g_k certain values if the conditions above
% are met or not the first element of the summation is if x is in a nodal
% region, the second element is if the x is in a myelinated region.
% a_fctn = @(x) a*a_fctn(x) + a_my*(a_fctn(x)==0); 
% C_m_fctn = @(x) C_m*C_m_fctn(x) + C_my*(C_m_fctn(x)==0);
% G_K_fctn = @(x) G_K_nodal*G_K_fctn(x) + G_K_internodal*(G_K_fctn(x)==0); 
% G_Na_fctn = @(x) G_Na_nodal*G_Na_fctn(x) + G_Na_internodal*(G_Na_fctn(x)==0); 


% Initialize Vm, n, m, and h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_0 = 0.6771; % probability that potassium gate is open 
M_0 = 0.4971; % probability that Sodium activation gate is open 
H_0 = 0.0512; % probability that Sodium inactivation gate is open
V_m0 = -40.1327; % (mV) Voltage 
Vm = V_m0 * ones(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);
N(1) = N_0;
M(1) = M_0;
H(1) = H_0;
for i = 2:m-1
    seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
    myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
    myelin_end = seg*(N_s); % End of internodal region in this segment

    if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
        newN(i) = 0;
        newM(i) = 0;
        newH(i) = 0;
    else % nodal region or end points 
        newN(i) = N_0;
        newM(i) = M_0;
        newH(i) = H_0;
    end
end



% defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

% Initialize matrix A and vectors b and f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(m);
g_1 = zeros(m, 1);
g_2 = zeros(m, 1);

% Starting the time loop
%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)
    
    % Updating the probability gate functions n, m and h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:m-1
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment

        if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
            newN(i) = 0;
            newM(i) = 0;
            newH(i) = 0;
        else % nodal region and end points 
            newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
            newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
            newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
        end
    end
    newN(m) = 0;
    newM(m) = 0;
    newH(m) = 0;

    % Constructing the A matrix
    A(1, 1) = 1;
    A(1, 2) = -1;
    A(m, m-1) = -1;
    A(m, m) = 1;
    % i is the space step
    for i = 2:m-1
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment
        seg_start = (seg - 1)*(N_s); % index of the start of the segment

        % Internodal region
        if (i > myelin_start + 1) && (i < myelin_end + 1)
            gamma1 = -dt*a_my/(2*R_i*C_my*dx^2);
            gamma2 = 1 + dt*a_my/(R_i*C_my*dx^2) + dt*G_K_internodal*(N(i)^4)/C_my + dt*G_Na_internodal*(M(i)^3)*H(i)/C_my + dt*G_L/C_my;
            gamma3 = -dt*a_my/(2*R_i*C_my*dx^2); 
            gamma4 = 1; 
            gamma5 = dt*G_K_internodal*N(i)^4*E_K/C_my + dt*G_Na_internodal*(M(i)^3)*H(i)*E_Na/C_my + dt*G_L*E_L/C_my;
        % nodal region and end points
        else
            gamma1 = -dt*a/(2*R_i*C_m*dx^2);
            gamma2 = 1 + dt*a/(R_i*C_m*dx^2) + dt*G_K_nodal*(N(i)^4)/C_m + dt*G_Na_nodal*(M(i)^3)*H(i)/C_m + dt*G_L/C_m;
            gamma3 = -dt*a/(2*R_i*C_m*dx^2); 
            gamma4 = 1; 
            gamma5 = dt*G_K_nodal*N(i)^4*E_K/C_m + dt*G_Na_nodal*(M(i)^3)*H(i)*E_Na/C_m + dt*G_L*E_L/C_m;
            % because the stimulus can only occur in the nodal region
            if (j*dt >= S_T0 && j*dt <= S_T1) && (i*dx >= S_P0 && i*dx <= S_P1)
                gamma2 = 1 + dt*a/(R_i*C_m*dx^2) + dt*G_K_nodal*(N(i)^4)/C_m + dt*(G_Na_nodal*(M(i)^3)*H(i) + S_v)/C_m + dt*G_L/C_m;
                gamma5 = dt*G_K_nodal*(N(i)^4)*E_K/C_m + dt*(G_Na_nodal*(M(i)^3)*H(i) + S_v)*G_Na_nodal/C_m + dt*G_L*E_L/C_m;
            end
        end

        % constructing the A and b matrix and vector
        A(i, i-1) = gamma1;
        A(i, i) = gamma2;
        A(i, i+1) = gamma3;
        g_1(i) = gamma4*Vm(i);
        g_2(i) = gamma5; 
       
    end

    % setting newU (the solution from Ax = b))
    newVm = transpose(A\(g_1+g_2));
    Vm = newVm;

    % adding the newly defined vectors to the 'all' matrices
    Vm_all(j+1,:) = Vm;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

    j

end


% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = L*0.25; % in cm
position2 = L*0.5; % in cm
position3 = L*0.75; % in cm
position4 = L; % in cm
position5 = 0.0002; 

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5];

% Times to observe the voltage along the axon
time1 = T*0.25; % in ms
time2 = T*0.5; % in ms
time3 = T*0.75; % in ms
time4 = T; % in ms


list_of_times = [time1
                 time2
                 time3
                 time4
                 ];

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

% save('HH_data.mat');
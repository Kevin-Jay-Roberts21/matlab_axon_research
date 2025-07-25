% Code used to simulate demylination and to construct nonuniform axon
% geometries
% Kevin Roberts
% July 2025

clear all
close all
clc

% Defining the constant Mesh Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.00005; % (cm) space step
dt = 0.01; % (ms) time step 
T = 30; % (ms) the total time of the experiment
n = T/dt + 1; % (#) n is the number of time steps

% Defining the constant material parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.55*10^(-4); % (cm) radius in nodal region
R_i = 0.0712; % (kilo-ohms*cm) intracellular resistivity
R_m = 24.8; % (kilo-ohms*cm^2) specific membrane resistance
C_m = 1.23; % (micro-farads/cm^2) specific membrane capacitance
G_K = 80; % (mS/cm^2) max specific potassium conductance
G_Na = 3000; % (mS/cm^2) max specific sodium conductance 
G_L = 80; % (mS/cm^2) specific leak conductance
E_K = -82; % (mV) Nernst potential for potassium ions
E_Na = 45; % (mV) Nernst potential for sodium ions
E_L = -59.4; % (mV) Nernst potential for leak channels
E_rest = -59.4; % (mV) effective resting nernst potential

% defining rho constant
rho = dt/dx^2;

% Stimulus Information
%%%%%%%%%%%%%%%%%%%%%%
S_v = 2000; % (mS/cm^2) % stimulus value
S_T0 = 1; % (ms) start time of when stimulus is added
S_T1 = 1.1; % (ms) end time of when stimulus is added 
S_P0 = 0; % (cm) start position of adding the stimulus (corresponds to ii = 1)
S_P1 = 0.0005; % (cm) end position of adding the stimulus (corresponds to ii = 11)
% in the S function ii, is the space index and tt is the time index (inclusive)
S = @(ii, tt) S_v * ((abs(tt * dt - S_T0) <= 1e-10 | tt * dt > S_T0) & ...
                    (tt * dt < S_T1 | abs(tt * dt - S_T1) <= 1e-10) & ...
                    (abs(ii * dx - S_P0) <= 1e-10 | ii * dx > S_P0) & ...
                    (ii * dx < S_P1 | abs(ii * dx - S_P1) <= 1e-10));

% Defining alpha/beta functions as well as the b_1, c_1 and f_1 functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_base = 20; % (C) base temperature
T_actual = 20; % (C) the temperature of the squid axon
Q_10_Na = 2.2; % (#) temperature coefficient for Na current
Q_10_K = 3; % (#) temperature coefficient for K current
phi_Na = Q_10_Na^((T_actual - T_base)/10); % (#) temperature scaling factor for Na current
phi_K = Q_10_K^((T_actual - T_base)/10); % (#) temperature scaling factor for K current
alpha_n = @(Vm) phi_K * 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10));
beta_n = @(Vm) phi_K * 0.125*exp(-(Vm + 65)/80);
alpha_m = @(Vm) phi_Na * 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10));
beta_m = @(Vm) phi_Na * 4*exp(-(Vm + 65)/18);
alpha_h = @(Vm) phi_Na * 0.07*exp(-(Vm + 65)/20);
beta_h = @(Vm) phi_Na * 1/(1 + exp(-(Vm + 35)/10));




% Handling other variables that may differ in each axon segment
% Note: Variables that will vary for each axon segment include L_n, L_my,
% a_my, R_my and C_my

% defining L_n, L_my, a_my, R_my and C_my in 20 axon segments
% Note: In this loop, a_my, R_my and C_my in axon segment 5 differs from the rest
for i = 1:20
    if i == 5
        segments(i).a_my = a/0.698*(0.9); % (cm) radius in myelinated region
        segments(i).R_my = 63.7/(0.9); % (kilo-ohms*cm^2) specfic myelin resistance
        segments(i).C_my = 0.113*(0.9); % (micro-fards/cm^2) specific myelin capacitance
        segments(i).L_n = 0.0005; % (cm) nodal length
        segments(i).L_my = 0.0075; % (cm) internodal length
    else
        segments(i).a_my = a/0.698; % (cm) radius in myelinated region
        segments(i).R_my = 63.7; % (kilo-ohms*cm^2) specfic myelin resistance
        segments(i).C_my = 0.113; % (micro-fards/cm^2) specific myelin capacitance
        segments(i).L_n = 0.0005; % (cm) nodal length
        segments(i).L_my = 0.0075; % (cm) internodal length
    end
end

n_s = length(segments); % number of axon segments

% Building the spatial grid
x = [];
region_type = [];
C_my_all = [];
R_my_all = [];
a_my_all = [];

for s = 1:n_s
    N_n = round(segments(s).L_n/dx);
    N_my = round(segments(s).L_my/dx);
    seg_length = N_n + N_my;

    for i = 1:seg_length
        x(end+1) = (length(x)) * dx;  % <--- FIXED: populate x

        if i == 1 || i == seg_length
            region_type(end+1) = 3; % endpoint
        elseif i <= N_n
            region_type(end+1) = 1; % nodal
        else
            region_type(end+1) = 2; % internodal
        end

        % Material property arrays
        C_my_all(end+1) = segments(s).C_my;
        R_my_all(end+1) = segments(s).R_my;
        a_my_all(end+1) = segments(s).a_my;
    end
end

m = length(x); % total number of space points

% Modifying b_1 and w_1 (no longer creating functions, defining vectors instead)
w_1_all = zeros(1, m);
b_1_all = zeros(1, m);

for i = 1:m
    if region_type(i) == 2  % Internodal
        w_1_all(i) = a^2/(C_my_all(i)*a_my_all(i)*R_i);
        b_1_all(i) = (a/(2*R_i*C_m))*(1 + (C_m*a)/(C_my_all(i)*a_my_all(i)));
    else                    % Nodal or end
        w_1_all(i) = 0;
        b_1_all(i) = a/(2*R_i*C_m);
    end
end

% Modifying f_1 and c_1 functions
% c_1 function
function val = c_1(i, n, m, h, tt, region_type, C_m, G_K, G_Na, G_L, S)
    if region_type(i) == 2  % internodal
        val = -1/(R_m*C_m);  % use R_m = 24.8 directly or pass as arg
    elseif region_type(i) == 1  % nodal
        val = -1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(i, tt)) + G_L);
    else  % endpoint (average)
        nodal_val = -1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(i, tt)) + G_L);
        internodal_val = -1/(R_m*C_m);
        val = (nodal_val + internodal_val)/2;
    end
end
 
% f_1 function
function val = f_1(i, Vmy, n, m, h, tt, region_type, C_my_all, R_my_all, E_rest, C_m, R_m, G_K, G_Na, G_L, E_K, E_Na, E_L, S)
    if region_type(i) == 2  % internodal
        term1 = (1/(R_m*C_m) - 1/(C_my_all(i)*R_my_all(i)))*Vmy;
        term2 = E_rest/(R_m*C_m);
        val = term1 + term2;
    elseif region_type(i) == 1  % nodal
        val = 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(i, tt))*E_Na + G_L*E_L);
    else  % endpoint (average of node + internode formula)
        internodal_term = (1/(R_m*C_m) - 1/(C_my_all(i)*R_my_all(i)))*Vmy + E_rest/(R_m*C_m);
        nodal_term = 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(i, tt))*E_Na + G_L*E_L);
        val = (internodal_term + nodal_term)/2;
    end
end


% Initialization
%%%%%%%%%%%%%%%%
V_m0 = -58.1539; % (mV) initial condition for membrane potential 
V_my0 = 0.8004; % (mV) initial condition for axon potential in periaxonal space
N_0 = 0.4258; % (dimless) initial condition for gating variable n
M_0 = 0.1144; % (dimless) initial condition for gating variable m
H_0 = 0.3560; % (dimless) initial condition for gating variable h
Vm = V_m0 * ones(1, m);
Vmy = zeros(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);
N(1) = N_0;
M(1) = M_0;
H(1) = H_0;
for i = 2:m-1
    if region_type(i) == 2 % internodal region
        Vmy(i) = V_my0;
        N(i) = 0;
        M(i) = 0;
        H(i) = 0;
    else % nodal and end points
        Vmy(i) = 0;
        N(i) = N_0;
        M(i) = M_0;
        H(i) = H_0;
    end
end

% Defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
Vm_minus_Vmy(1,:) = Vm - Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

% Running the time loop
%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)
    
    % Updating the probability gate functions n, m and h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newN = zeros(1, m);
    newM = zeros(1, m);
    newH = zeros(1, m);

    for i = 2:(m-1)
        if region_type(i) == 2 % internodal
            newN(i) = 0;
            newM(i) = 0;
            newH(i) = 0;
        else % nodal or endpoint
            newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
            newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
            newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
        end
    end
    newN(m) = 0;
    newM(m) = 0;
    newH(m) = 0;
    
    % Updating Vmy
    %%%%%%%%%%%%%%
    newVmy = zeros(1, m);
    newVmy(1) = 0;

    for i = 2:(m-1)
        if region_type(i) == 2  % internodal region
            eta1 = 1 / (1 + dt / (C_my_all(i) * R_my_all(i)));
            eta2 = rho * w_1_all(i) / 2 * eta1;
            eta3 = -rho * w_1_all(i) * eta1;
            eta4 = rho * w_1_all(i) / 2 * eta1;
        else % nodal and endpoint regions
            eta1 = 0;
            eta2 = 0;
            eta3 = 0;
            eta4 = 0;
        end
        newVmy(i) = eta1*Vmy(i) + eta2*Vm(i-1) + eta3*Vm(i) + eta4*Vm(i+1);
    end
    newVmy(m) = 0;
    
    % Updating Vm
    %%%%%%%%%%%%%
    % Defining the A matrix and RHS vectors
    A = zeros(m, m);
    g_1 = zeros(m, 1);
    g_2 = zeros(m, 1);

    % Apply Neumann BC at boundaries (dV/dx = 0)
    A(1, 1) = 1;
    A(1, 2) = -1;
    A(m, m-1) = -1;
    A(m, m) = 1;

    for i = 2:(m-1)
        % Approximate b(i ± 1/2) using central difference
        b_i_minus = (b_1_all(i-1) + b_1_all(i)) / 2;
        b_i_plus  = (b_1_all(i+1) + b_1_all(i)) / 2;

        gamma1 = -rho * b_i_minus;
        gamma2 = 1 - dt * c_1(i, newN(i), newM(i), newH(i), j, region_type, C_m, G_K, G_Na, G_L, S) ...
                     + rho * (b_i_minus + b_i_plus);
        gamma3 = -rho * b_i_plus;
        gamma4 = 1;
        gamma5 = dt * f_1(i, newVmy(i), newN(i), newM(i), newH(i), j, ...
                         region_type, C_my_all, R_my_all, E_rest, ...
                         C_m, R_m, G_K, G_Na, G_L, E_K, E_Na, E_L, S);

        A(i, i-1) = gamma1;
        A(i, i)   = gamma2;
        A(i, i+1) = gamma3;
        g_1(i) = gamma4 * Vm(i);
        g_2(i) = gamma5;
    end

    % Solve A * Vm^{j+1} = g_1 + g_2
    newVm = (A \ (g_1 + g_2))';

    % Update all variables
    Vm = newVm;
    Vmy = newVmy;
    N = newN;
    M = newM;
    H = newH;

    % Store data
    Vm_all(j+1,:) = Vm;
    Vmy_all(j+1,:) = Vmy;
    Vm_minus_Vmy(j+1,:) = Vm - Vmy;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

    % Progress update
    if mod(j, 50) == 0
        fprintf("Completed time step %d / %d\n", j, n-1);
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PICKING TIME AND POSITION SHOTS TO PLOT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING Vm SPATIAL AND TEMPORAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First figure: Voltage along the axon at different times
figure(1);

% Adjust the figure size (Position [left, bottom, width, height])
set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)

% Create subplot (1 row, 2 columns, 1st subplot)
subplot(1, 2, 1);
t1 = linspace(0, 10000*L, m);
plot(t1, Vm_all(round(time1/dt),:));
hold on
for i = 2:length(list_of_times)
    plot(t1, Vm_all(round(list_of_times(i)/dt),:));
end

% Describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('$V_m$ at t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex');
ylabel("$V_m$ in millivolts.", 'Interpreter','latex');
xlabel("Length of the axon in um.");


% Second figure: Voltage vs Time at different positions
% Create subplot (1 row, 2 columns, 2nd subplot)
subplot(1, 2, 2);
t2 = linspace(0, T, n); % FULL MATRIX
plot(t2, Vm_all(:,round(position1/dx)));
hold on
for i = 2:length(list_of_positions)
    plot(t2, Vm_all(:,round(list_of_positions(i)/dx)));
end

% Describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('$V_m$ at x = %g um', 10000*list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex');
ylabel("$V_m$ in millivolts.", 'Interpreter', 'latex');
xlabel("Time in milliseconds.");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING Vmy SPATIAL AND TEMPORAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First figure: Voltage along the axon at different times
figure(2);

% Adjust the figure size (Position [left, bottom, width, height])
set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)

% Create subplot (1 row, 2 columns, 1st subplot)
subplot(1, 2, 1);
t1 = linspace(0, L*10000, m);
plot(t1, Vmy_all(round(time1/dt),:));
hold on
for i = 2:length(list_of_times)
    plot(t1, Vmy_all(round(list_of_times(i)/dt),:));
end

% Describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('$V_{my}$ at t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex');
ylabel("$V_{my}$ in millivolts.", 'Interpreter','latex');
xlabel("Length of the axon in um.");


% Second figure: Voltage vs Time at different positions
% Create subplot (1 row, 2 columns, 2nd subplot)
subplot(1, 2, 2);
t2 = linspace(0, T, n); % FULL MATRIX
plot(t2, Vmy_all(:,round(position1/dx)));
hold on
for i = 2:length(list_of_positions)
    plot(t2, Vmy_all(:,round(list_of_positions(i)/dx)));
end

% Describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('$V_{my}$ at x = %g um', 10000*list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex');
ylabel("$V_{my}$ in millivolts.", 'Interpreter', 'latex');
xlabel("Time in milliseconds.");




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING Vm-Vmy SPATIAL AND TEMPORAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First figure: Voltage along the axon at different times
figure(3);

% Adjust the figure size (Position [left, bottom, width, height])
set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)

% Create subplot (1 row, 2 columns, 1st subplot)
subplot(1, 2, 1);
t1 = linspace(0, 10000*L, m);
plot(t1, Vm_minus_Vmy(round(time1/dt),:));
hold on
for i = 2:length(list_of_times)
    plot(t1, Vm_minus_Vmy(round(list_of_times(i)/dt),:));
end

% Describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('$V_m - V_{my}$ at t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex');
ylabel("$V_m - V_{my}$ in millivolts.", 'Interpreter','latex');
xlabel("Length of the axon in um.");


% Second figure: Voltage vs Time at different positions
% Create subplot (1 row, 2 columns, 2nd subplot)
subplot(1, 2, 2);
t2 = linspace(0, T, n); % FULL MATRIX
plot(t2, Vm_minus_Vmy(:,round(position1/dx)));
hold on
for i = 2:length(list_of_positions)
    plot(t2, Vm_minus_Vmy(:,round(list_of_positions(i)/dx)));
end

% Describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('$V_m - V_{my}$ at x = %g um', 10000*list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex');
ylabel("$V_m - V_{my}$ in millivolts.", 'Interpreter', 'latex');
xlabel("Time in milliseconds.");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING n, m, and h TEMPROAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
plot(t2, N_all(:,round(position5/dx)))
hold on
plot(t2, M_all(:,round(position5/dx)))
hold on
plot(t2, H_all(:,round(position5/dx)))
legendStrings3 = {
    sprintf('n at x = %g um', 10000*position5), ...
    sprintf('m at x = %g um', 10000*position5), ...
    sprintf('h at x = %g um', 10000*position5)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")

% save('SC_Cohen_set1_T66.mat'); 

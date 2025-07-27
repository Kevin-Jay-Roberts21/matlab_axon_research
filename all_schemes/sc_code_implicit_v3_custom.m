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
T = 10; % (ms) the total time of the experiment
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

% Axon Segments
%%%%%%%%%%%%%%%
for i = 1:20
    if i == 5
        a_my = a/0.698; % a/0.698*0.9; 
        R_my = 63.7; % 63.7/0.9; 
        C_my = 0.113; % 0.113*0.9;
    else
        a_my = a/0.698; 
        R_my = 63.7; 
        C_my = 0.113;
    end
    segments(i) = struct('a_my', a_my, 'R_my', R_my, 'C_my', C_my, 'L_n', 0.0005, 'L_my', 0.0075);
end

% Spatial Grid
%%%%%%%%%%%%%%
x = [];
region_type = [];
C_my_split = [];
R_my_split = [];
a_my_split = [];

% Add initial endpoint at the very start of the axon
x(end+1) = 0;
region_type(end+1) = 3;  % global starting endpoint
C_my_split(end+1) = segments(1).C_my;
R_my_split(end+1) = segments(1).R_my;
a_my_split(end+1) = segments(1).a_my;

for seg = 1:length(segments)
    N_n = round(segments(seg).L_n / dx);
    N_my = round(segments(seg).L_my / dx);

    % Nodal region
    for i = 1:(N_n-1)
        x(end+1) = length(x) * dx;
        region_type(end+1) = 1;
        C_my_split(end+1) = segments(seg).C_my;
        R_my_split(end+1) = segments(seg).R_my;
        a_my_split(end+1) = segments(seg).a_my;
    end

    % Endpoint between nodal and internodal
    x(end+1) = length(x) * dx;
    region_type(end+1) = 3;
    C_my_split(end+1) = segments(seg).C_my;
    R_my_split(end+1) = segments(seg).R_my;
    a_my_split(end+1) = segments(seg).a_my;

    % Internodal region
    for i = 1:(N_my-1)
        x(end+1) = length(x) * dx;
        region_type(end+1) = 2;
        C_my_split(end+1) = segments(seg).C_my;
        R_my_split(end+1) = segments(seg).R_my;
        a_my_split(end+1) = segments(seg).a_my;
    end

    % Endpoint at end of segment
    x(end+1) = length(x) * dx;
    region_type(end+1) = 3;
    C_my_split(end+1) = segments(seg).C_my;
    R_my_split(end+1) = segments(seg).R_my;
    a_my_split(end+1) = segments(seg).a_my;
end

m = length(x);
L = (m)*dx;

% Coefficient Arrays
%%%%%%%%%%%%%%%%%%%%
w_1_split = zeros(1, m);
b_1_split = zeros(1, m);
for i = 1:m
    if region_type(i) == 2
        w_1_split(i) = a^2/(C_my_split(i)*a_my_split(i)*R_i);
        b_1_split(i) = (a/(2*R_i*C_m)) * (1 + (C_m*a)/(C_my_split(i)*a_my_split(i)));
    else
        w_1_split(i) = 0;
        b_1_split(i) = a/(2*R_i*C_m);
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
for i = 1:m
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
Vm_minus_Vmy(1,:) = Vm - Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;


% Defining c_1 and f_1 functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_1 function
c_1 = @(ii, n, m, h, tt)...
    (region_type(i) == 2)*(-1/(R_m*C_m)) + ... % internodal
    (region_type(i) == 1)*(-1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(ii, tt)) + G_L)) + ... % nodal
    (region_type(i) == 3)*((-1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(ii, tt)) + G_L) + -1/(R_m*C_m))/2); % end point
 
% f_1 function
f_1 = @(ii, Vmy, n, m, h, tt)...
    (region_type(i) == 2)*((1/(R_m*C_m) - 1/(C_my_split(i)*R_my_split(i)))*Vmy + E_rest/(R_m*C_m)) + ... % internodal
    (region_type(i) == 1)*(1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(i, tt))*E_Na + G_L*E_L)) + ... % nodal
    (region_type(i) == 3)*(((1/(R_m*C_m) - 1/(C_my_split(i)*R_my_split(i)))*Vmy + E_rest/(R_m*C_m) + 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(i, tt))*E_Na + G_L*E_L))/2); % end point


% Running the time loop
%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)
    % Update n, m, h
    for i = 1:m
        type = region_type(i);
        switch type
            case 1  % nodal
                newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
                newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
                newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
            case 2  % internodal
                newN(i) = 0;
                newM(i) = 0;
                newH(i) = 0;
            otherwise  % endpoint
                newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
                newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
                newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
        end
    end

    % Update Vmy
    newVmy(1) = 0;
    for i = 2:m-1
        if region_type(i) == 2  % internodal
            w1 = w_1_split(i);
            eta1 = 1/(1 + dt/(C_my_split(i)*R_my_split(i)));
            eta2 = rho * w1/2 * eta1;
            eta3 = -rho * w1 * eta1;
            eta4 = rho * w1/2 * eta1;
        else
            eta1 = 0;
            eta2 = 0;
            eta3 = 0;
            eta4 = 0;
        end
        newVmy(i) = eta1*Vmy(i) + eta2*Vm(i-1) + eta3*Vm(i) + eta4*Vm(i+1);
    end
    newVmy(m) = 0;

    % Assemble A matrix, g_1, g_2
    A = zeros(m);
    g_1 = zeros(m,1);
    g_2 = zeros(m,1);
    A(1,1) = 1; A(1,2) = -1;
    A(m,m-1) = -1; A(m,m) = 1;

    for i = 2:m-1
        bL = b_1_split(i-1);
        bC = b_1_split(i);
        bR = b_1_split(i+1);
        gamma1 = -rho*bL;
        gamma2 = 1 - dt*c_1(i, newN(i), newM(i), newH(i), j) + rho*(bL + bR);
        gamma3 = -rho*bR;
        gamma4 = 1;
        gamma5 = dt * f_1(i, newVmy(i), newN(i), newM(i), newH(i), j);

        A(i,i-1) = gamma1;
        A(i,i) = gamma2;
        A(i,i+1) = gamma3;
        g_1(i) = gamma4*Vm(i);
        g_2(i) = gamma5;
    end

    newVm = transpose(A\(g_1 + g_2));

    Vm = newVm;
    Vmy = newVmy;
    N = newN;
    M = newM;
    H = newH;

    Vm_all(j+1,:) = Vm;
    Vmy_all(j+1,:) = Vmy;
    Vm_minus_Vmy(j+1,:) = Vm - Vmy;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

    j
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

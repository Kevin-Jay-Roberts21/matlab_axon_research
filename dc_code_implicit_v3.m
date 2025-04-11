% Double Cable Model Finite Difference implicit discretization
% Kevin Roberts
% February 2025

clear all
close all
clc

% Defining the Thickness, Length and other Mesh Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.00005; % (cm) space step
dt = 0.01; % (ms) time step 
L_my = 0.0075; % (cm) internodal length
L_n = 0.0005; % (cm) nodal length
L_pn = 2.3*10^(-4); % (cm) paranodal length
d_pa = 12.3*10^(-7); % (cm) periaxonal thickness
d_pn = 7.4*10^(-7); % (cm) paranodal thickness
L_s = L_n + L_my; % (cm) length of an axon segment
n_s = 20; % (dimless) number of axon segments
L = n_s*L_s; % (cm) total length of axon
T = 30; % (ms) the total time of the experiment
N_n = round(L_n/dx); % number of space steps in a nodal region
N_my = round(L_my/dx); % number of space steps in an internodal region
N_s = N_n + N_my; % number of space steps in an entire axon segement
m = N_s*n_s + 1; % total number of space steps
n = T/dt + 1; % n is the number of time steps

% Defining the material properties on other intrinsic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to change to test Dr. Huang's Paper
% a_my = 0.0001; % (cm) radius in myelinated region
% a = 0.0001238; % (cm) radius in nodal region
% R_my = 123.6795; % (kilo-ohms*cm^2) specfic myelin resistance
% C_my = 0.0081; % (micro-fards/cm^2) specific myelin capacitance

% original a, a_my, R_my and C_my
a = 0.55*10^(-4); % (cm) axon radius in nodal region
a_my = a/0.698; % (cm) axon radius in myelinated section 
R_i = 0.155; % (kilo-ohms*cm) intracellular resistivity
R_m = 24.6; % (kilo-ohms*cm^2) specific membrane resistance
C_m = 1.15; % (micro-farads/cm^2) specific membrane capacitance
r_pa = 125*10^6; % 96.3*10^6; % (kilo-ohms/cm) periaxonal resistivity per unit length
R_pa = r_pa*pi*d_pa*(2*a + d_pa); % (kilo-ohms*cm) resistivity of the periaxonal space (computed)
r_pn = 2450*10^6; % (kilo-ohms/cm) paranodal resitance per unit length (used in BC since r_bar_pn = r_pn * L_pn) 
R_my = 240; % (kilo-ohms*cm^2) specific myelin resistance
C_my = 0.0379; % (micro-farads/cm^2) specific myelin capacitance
G_K = 80; % (mS/cm^2) max specific potassium conductance
G_Na = 3000; % (mS/cm^2) max specific sodium conductance 
G_L = 80; % (mS/cm^2) specific leak conductance
E_K = -82; % (mV) Nernst potential for potassium ions
E_Na = 45; % (mV) Nernst potential for sodium ions
E_L = -59.4; % (mV) Nernst potential for leak channels
E_rest = -59.4; % (mV) effective resting nernst potential

% defining rho, w1, w2 and w3 constants
rho = dt/dx^2; % creating the courant number
w1 = a^2/(C_my*a_my*R_i);
w2 = d_pa*(2*a + d_pa)/(C_my*a_my*R_pa);
w3 = r_pa/(r_pn*L_pn);

% Stimulus Information
%%%%%%%%%%%%%%%%%%%%%%
S_v = 970; % (in mS/cm^2) % stimulus value
S_T0 = 1; % start time of when stimulus is added (in ms)
S_T1 = 1.1; % end time of when stimulus is added (in ms)
S_P0 = 0.0001; % start position of adding the stimulus (in cm)
S_P1 = 0.0004; % end position of adding the stimulus (in cm)
% in the S functi on ii, is the space index and tt is the time index
S = @(ii, tt) S_v * ((abs(tt * dt - S_T0) <= 1e-10 | tt * dt > S_T0) & ...
                    (tt * dt < S_T1 | abs(tt * dt - S_T1) <= 1e-10) & ...
                    (abs(ii * dx - S_P0) <= 1e-10 | ii * dx > S_P0) & ...
                    (ii * dx < S_P1 | abs(ii * dx - S_P1) <= 1e-10));

% Defining alpha/beta functions as well as the b_1, c_1 and f_1 functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_base = 20; % (C) base temperature
T_actual = 20; % (C) the temperature of the squid axon
Q_10_Na = 2.2; % (dimless) temperature coefficient for Na current
Q_10_K = 3; % (dimless) temperature coefficient for K current
phi_Na = Q_10_Na^((T_actual - T_base)/10); % (dimless) temperature scaling factor for Na current
phi_K = Q_10_K^((T_actual - T_base)/10); % (dimless) temperature scaling factor for K current
alpha_n = @(Vm) phi_K * 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10));
beta_n = @(Vm) phi_K * 0.125*exp(-(Vm + 65)/80);
alpha_m = @(Vm) phi_Na * 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10));
beta_m = @(Vm) phi_Na * 4*exp(-(Vm + 65)/18);
alpha_h = @(Vm) phi_Na * 0.07*exp(-(Vm + 65)/20);
beta_h = @(Vm) phi_Na * 1/(1 + exp(-(Vm + 35)/10));

% defining the b_1(x_i) function
B_1 = (a/(2*R_i*C_m))*(1 + C_m*a/(C_my*a_my)); % Internodal region
B_2 = a/(2*R_i*C_m); % Nodal region
B_3 = (B_1 + B_2)/2; % End point (may not be used)
b_1 = @(ii) (mod(ii - 1, N_s) > N_n).*B_1 + ... % Internodal region
           (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*B_2 + ... % Nodal region
           ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*B_3; % End point      

% defining the c_1(x_i) function
C_1 = -1/(R_m*C_m); % Internodal region
C_2 = @(n, m, h, ii, tt) -1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(ii, tt)) + G_L); % Nodal region
C_3 = @(n, m, h, ii, tt) (C_1 + C_2(n, m, h, ii, tt))/2; % End point
c_1 = @(n, m, h, ii, tt) (mod(ii - 1, N_s) > N_n).*C_1 + ... % Internodal region
           (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*C_2(n, m, h, ii, tt) + ... % Nodal region
           ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*C_3(n, m, h, ii, tt); % End point        
       
% defining the f_1(x_i) function
F_1 = @(Vmy_i) (1/(R_m*C_m) - 1/(C_my*R_my))*Vmy_i + E_rest/(R_m*C_m);
F_4 = @(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1) w2/(2*dx^2)*Vmy_i_minus_1 - w2/dx^2*Vmy_i + w2/(2*dx^2)*Vmy_i_plus_1; % Internodal region
F_2 = @(n, m, h, ii, tt) 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(ii, tt))*E_Na + G_L*E_L); % Nodal region
F_5 = @(Vmy_i, n, m, h, ii, tt) (F_1(Vmy_i) + F_2(n, m, h, ii, tt))/2; % End point
 
f_2 = @(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1, n, m, h, ii, tt) (mod(ii - 1, N_s) > N_n).*(F_1(Vmy_i) + F_4(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1)) + ... % Internodal region
           (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*F_2(n, m, h, ii, tt) + ... % Nodal region
           ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*F_5(Vmy_i, n, m, h, ii, tt); % End point

% Initialization
%%%%%%%%%%%%%%%%
V_m0 = -58.1124; % (mV) initial condition for membrane potential 
V_my0 = 0.03; % (mV) initial condition for axon potential in periaxonal space
N_0 = 0.4264; % (dimless) initial condition for gating variable n
M_0 = 0.1148; % (dimless) initial condition for gating variable m
H_0 = 0.3548; % (dimless) initial condition for gating variable h
Vm = V_m0 * ones(1, m);
Vmy = V_my0 * ones(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);
Vmy(1) = 0;
N(1) = N_0;
M(1) = M_0;
H(1) = H_0;
for i = 2:m-1 % because we want to keep the end points (1 and M) for Vmy, n, m, h at 0
    
    seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
    myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
    myelin_end = seg*(N_s); % End of internodal region in this segment
    seg_start = (seg - 1)*(N_s); % index of the start of the segment

    % Internodal region
    if (i > myelin_start + 1) && (i < myelin_end + 1)
        Vmy(i) = V_my0;
        N(i) = 0;
        M(i) = 0;
        H(i) = 0;
    % Nodal region
    elseif (i > seg_start + 1) && (i < myelin_start + 1)
        Vmy(i) = 0;
        N(i) = N_0;
        M(i) = M_0; 
        H(i) = H_0;
    % End points
    elseif (i == seg_start + 1) % Right end point (x_R)
        Vmy(i) = V_my0/(1 + w3*dx);
        N(i) = N_0;
        M(i) = M_0;
        H(i) = H_0;
    elseif (i == myelin_start + 1) % Left end point (x_L)
        Vmy(i) = V_my0/(1 + w3*dx);
        N(i) = N_0;
        M(i) = M_0;
        H(i) = H_0;
    end 
end
Vmy(m) = V_my0/(1 + w3*dx); % setting the right end point of the last myelinated section of the axon

% defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
Vm_minus_Vmy(1,:) = Vm - Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

% Running the time loop
%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)
    
    % updating the probability gate functions n, m and h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:m-1
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment
        seg_start = (seg - 1)*(N_s); % index of the start of the segment

        if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
            newN(i) = 0;
            newM(i) = 0;
            newH(i) = 0;
        elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
            newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
            newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
            newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
        else % End point
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
    
    % Defining the A_2 matrix as well as e_2 and f_2
    A_2 = zeros(m, m);
    e_1 = zeros(m, 1);
    e_2 = zeros(m, 1);

    % using the boundary conditions to define the top and bottom row of A_2
    A_2(1, 1) = 1;
    A_2(1, 2) = 0;
    A_2(m, m-1) = -1;
    A_2(m, m) = (1 + w3*dx);

    for i = 2:(m-1)
        % updating the A_2 matrix each row is determined based on
        % internodal, nodal, and end points (left or right end points)
        seg = floor((i - 1)/N_s) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*N_s + N_n; % Start of internodal region in this segment
        myelin_end = seg*N_s; % End of internodal region in this segment
        seg_start = (seg - 1)*N_s; % index of the start of the segment

        if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
            eta1 = -rho*w2/2;
            eta2 = 1 + dt/(R_my*C_my) + rho*w2;
            eta3 = -rho*w2/2; 
            eta4 = 1; 
            eta5 = rho*w1/2*Vm(i-1) - rho*w1*Vm(i) + rho*w1/2*Vm(i+1);
        elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
            eta1 = 0;
            eta2 = 1;
            eta3 = 0; 
            eta4 = 0; 
            eta5 = 0;
        elseif (i == seg_start + 1) % Right end point (x_R)
            eta1 = -1;
            eta2 = (1 + dx*w3);
            eta3 = 0; 
            eta4 = 0; 
            eta5 = 0;
        elseif (i == myelin_start + 1) % Left end point (x_L)
            eta1 = 0;
            eta2 = (1 + dx*w3);
            eta3 = -1;
            eta4 = 0; 
            eta5 = 0;
        end

        A_2(i, i-1) = eta1;
        A_2(i, i) = eta2;
        A_2(i, i+1) = eta3;
        e_1(i) = eta4*Vmy(i);
        e_2(i) = eta5;

    end
    newVmy = transpose(A_2\(e_1+e_2));
    
    % Updating Vm 
    %%%%%%%%%%%%%

    % Defining the A_1 matrix
    A_1 = zeros(m, m);
    g_1 = zeros(m, 1);
    g_2 = zeros(m, 1);

    % using the boundary conditions to define the top and bottom row of A
    A_1(1, 1) = 1;
    A_1(1, 2) = -1;
    A_1(m, m-1) = -1;
    A_1(m, m) = 1;
    
    
    for i = 2:(m-1)

        gamma1 = -rho*b_1(i - 1/2);
        gamma2 = 1 - dt*c_1(newN(i), newM(i), newH(i), i, j) + rho*(b_1(i + 1/2) + b_1(i - 1/2));
        gamma3 = -rho*b_1(i + 1/2);
        gamma4 = 1;
        gamma5 = dt * f_2(newVmy(i-1), newVmy(i), newVmy(i+1), newN(i), newM(i), newH(i), i, j);

        A_1(i, i-1) = gamma1;
        A_1(i, i) = gamma2;
        A_1(i, i+1) = gamma3;
        g_1(i) = gamma4*Vm(i);
        g_2(i) = gamma5;
        
    end
    
    j % showing the j index, just for seeing how long simulation takes     

    % Solving for V_m^{j+1}
    %%%%%%%%%%%%%%%%%%%%%%%
    newVm = transpose(A_1\(g_1+g_2));

    % updating Vmy and Vm and adding the data to the _all matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vm = newVm;
    Vmy = newVmy;
    N = newN;
    M = newM;
    H = newH;
    
    % Adding the updated vectors to the 'all' data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vm_all(j+1,:) = Vm;
    Vmy_all(j+1,:) = Vmy;
    Vm_minus_Vmy(j+1,:) = Vm - Vmy;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

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

% save('DC_Cohen_DC_avg_params.mat');
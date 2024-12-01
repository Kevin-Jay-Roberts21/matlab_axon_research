% Solving the Single Cable Model Using Finite Difference Method
% Kevin Roberts
% November 2024

% SPECIAL NOTE: many of the for loops do the same thing, it may be worth
% combining some of these for loops later

clear all 
close all
clc

% Defining the material properties on other intrinsic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 5.5*10^(-5); % (cm) axon radius in nodal region
a_my = 5.623*10^(-5); % (cm) axon radius in myelinated section 
C_m = 0.001; % (ms/(ohms*cm^2)) specific membrane capacitance
C_my = 1.66*10^(-4); % (ms/(ohms*cm^2)) specific myelin capacitance
R_i = 150; % (ohms*cm) intracellular resistivity
R_my = 2.4*10^5; % (ohms*cm^2) specific myelin resistance
R_m = 2.5*10^3; % (ohms*cm^2) specific membrane resistance
G_K = 0.036; % (S/cm^2) max specific potassium conductance
G_Na = 0.12; % (S/cm^2) max specific sodium conductance 
G_L = 0.0003; % (S/cm^2) specific leak conductance
E_K = -82; % (mV) Nernst potential for potassium ions
E_Na = 45; % (mV) Nernst potential for sodium ions
E_L = -59.4; % (mV) Nernst potential for leak channels

% Defining the Mesh Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.0001; % (cm) space step
dt = 0.01; % (ms) time step 
L_my = 0.0075; % (cm) internodal length
L_n = 0.0005; % (cm) nodal length
L_s = L_n + L_my; % (cm) length of an axon segment
n_s = 10; % (dimless) number of axon segments
L = n_s*L_s; % (cm) total length of axon
T = 10; % (ms) the total time of the experiment
N_n = round(L_n/dx); % number of space steps in a nodal region
N_my = round(L_my/dx); % number of space steps in an internodal region
N_s = N_n + N_my; % number of space steps in an entire axon segement
m = L/dx + 1; % total number of space steps
n = T/dt; % n is the number of time steps

% Defining alpha/beta functions as well as the b_1 and f_1 functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defining the alpha_xi and beta_xi functions
alpha_n = @(Vm) 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10));
beta_n = @(Vm) 0.125*exp(-(Vm + 65)/80);
alpha_m = @(Vm) 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10));
beta_m = @(Vm) 4*exp(-(Vm + 65)/18);
alpha_h = @(Vm) 0.07*exp(-(Vm + 65)/20);
beta_h = @(Vm) 1/(1 + exp(-(Vm + 35)/10));

% defining the b_1(x_i)
B_1 = (a/(2*R_i))*(1 + C_m*a/(C_my*a_my)); % for internodal
B_2 = a/(2*R_i); % for nodal
B_3 = (B_1 + B_2)/2; % for end point
b_1 = @(ii) (mod(ii - 1, N_s) > N_n).*B_1 + ... % Internodal region
           (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*B_2 + ... % Nodal region
           ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*B_3; % Boundary point      

% defining the f_1(x_i) function
F_1 = @(Vm, Vmy) -Vm/R_m + (1/R_m - C_m/(C_my*R_my))*Vmy; % for internodal
F_2 = @(Vm, n, m, h) (G_K*n^4 - G_Na*m^3*h - G_L)*Vm + G_K*n^4*E_L + G_Na*m^3*h*E_Na + G_L*E_L; % for nodal
F_3 = @(Vm, Vmy, n, m, h) (F_1(Vm, Vmy) + F_2(Vm, n, m, h))/2; % for end point
f_1 = @(ii, Vm, Vmy, n, m, h) (mod(ii - 1, N_s) > N_n).*F_1(Vm, Vmy) + ... % Internodal region
           (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*F_2(Vm, n, m, h) + ... % Nodal region
           ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*F_3(Vm, Vmy, n, m, h); % Boundary point        


% Initialization
%%%%%%%%%%%%%%%%
V_m0 = -64.999; % (mV) initial condition for membrane potential 
V_my0 = -1; % (mV) initial condition for axon potential in periaxonal space
N_0 = 0.3177; % (dimless) initial condition for gating variable n
M_0 = 0.0529; % (dimless) initial condition for gating variable m
H_0 = 0.5961; % (dimless) initial condition for gating variable h
Vm = V_m0 * ones(1, m);
Vmy = zeros(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);
for i = 2:m-1 % because we want to keep the end points for Vmy, n, m, h at 0
    % regions of N, M and H
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
    % End points (we let Vmy, and N,M,H have ne effect on the end points)
    else
        Vmy(i) = V_my0;
        N(i) = N_0;
        M(i) = M_0;
        H(i) = H_0;
    end 
end

% defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
Vm_minus_Vmy(1,:) = Vm - Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;


% Defining the A matrix
%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(m, m);
% using the boundary conditions to define the top and bottom row of A
A(1, 1) = 1/dx;
A(1, 2) = -1/dx;
A(m, m-1) = -1/dx;
A(m, m) = 1/dx;

for i = 2:(m-1)
    % updating the A matrix rows if in an internodal region
    % NOTE: the b function accounts for the piecewise values
    % putting 0's into the nodal regions of Vmy and 0's into the internodal
    % regions of N, M and H
    seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
    myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
    myelin_end = seg*(N_s); % End of internodal region in this segment
    seg_start = (seg - 1)*(N_s); % index of the start of the segment

    % If i is in an Internodal Region
    if (i > myelin_start + 1) && (i < myelin_end + 1)
        A(i, i-1) = -b_1(i)/dx^2;
        A(i, i) = C_m/dt + 2*b_1(i)/dx^2;
        A(i, i+1) = -b_1(i)/dx^2;
    % If i is in a Nodal Region
    elseif (i > seg_start + 1) && (i < myelin_start + 1)
        A(i, i-1) = -b_1(i)/dx^2;
        A(i, i) = C_m/dt + 2*b_1(i)/dx^2;
        A(i, i+1) = -b_1(i)/dx^2;
    % If i is at an End Point (MIGHT HAVE TO SPECIFY FOR LEFT OR RIGHT END POINT)
    else 
        A(i, i-1) = -B_1/dx^2;
        A(i, i) = C_m/dt + 2*B_3/dx^2;
        A(i, i+1) = -B_2/dx^2;
    end
end

% Running the time loop
%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:(n-1)
    
    % updating Vmy
    for i = 2:m-1
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment
        seg_start = (seg - 1)*(N_s); % index of the start of the segment

        if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
            newVmy(i) = dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i-1) - dt*a^2/(C_my*a_my*R_i*dx^2)*Vm(i) + dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i+1) + (C_my/dt - 1/R_my)*Vmy(i);
        elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
            newVmy(i) = 0;
        else % End point
            newVmy(i) = dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i-1) - dt*a^2/(C_my*a_my*R_i*dx^2)*Vm(i) + dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i+1) + (C_my/dt - 1/R_my)*Vmy(i);
        end
    end
    
    % updating the probability gate functions n, m and h
    for i = 2:m-1
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment
        seg_start = (seg - 1)*(N_s); % index of the start of the segment

        if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
            newN(i) = 0;
            newM(i) = 0;
            newH(i) = 0;
        elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
            newN(i) = 1/(1/dt + alpha_n(Vm(i)) + beta_n(Vm(i))) * (N(i)/dt + alpha_n(Vm(i)));
            newM(i) = 1/(1/dt + alpha_m(Vm(i)) + beta_m(Vm(i))) * (M(i)/dt + alpha_m(Vm(i)));
            newH(i) = 1/(1/dt + alpha_h(Vm(i)) + beta_h(Vm(i))) * (H(i)/dt + alpha_h(Vm(i)));
        else % End point
            newN(i) = 1/(1/dt + alpha_n(Vm(i)) + beta_n(Vm(i))) * (N(i)/dt + alpha_n(Vm(i)));
            newM(i) = 1/(1/dt + alpha_m(Vm(i)) + beta_m(Vm(i))) * (M(i)/dt + alpha_m(Vm(i)));
            newH(i) = 1/(1/dt + alpha_h(Vm(i)) + beta_h(Vm(i))) * (H(i)/dt + alpha_h(Vm(i)));
        end
    end
    
    % updating the f function (end points are 0)
    f = zeros(m, 1);
    for i = 2:m-1
        f(i, 1) = f_1(i, Vm(i), Vmy(i), N(i), M(i), H(i)); 
    end
    % updating the b function (end points are 0)
    b = zeros(m, 1);
    for i = 2:m-1
        b(i, 1) = C_m/dt*Vm(i); 
    end
    
    j % showing the j index, just for seeing how long simulation takes 

    % because this is technically in a nodal region at the end point
    % recall this experiment is not symmetric
    newVmy(m) = 0;
    newN(m) = 0;
    newM(m) = 0;
    newH(m) = 0;

    % Solving for V_m^{j+1}
    newVm = transpose(A\(b+f));

    % updating Vmy and Vm and adding the data to the _all matrices
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

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PICKING TIME AND POSITION SHOTS TO PLOT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

position1 = 0.0245; % in cm
position2 = 0.03; % in cm
position3 = 0.0325; % in cm
position4 = 0.05; % in cm
position5 = 0.06; % in cm
position6 = 0.07; % in cm
position7 = 0.08; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     ];

% Times to observe the voltage along the axon
time1 = 1; % in ms
time2 = 2; % in ms
time3 = 4; % in ms
time4 = 5; % in ms
time5 = 7; % in ms
time6 = 8; % in ms
time7 = 10; % in ms

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
plot(t2, N_all(:,round(position3/dx)))
hold on
plot(t2, M_all(:,round(position3/dx)))
hold on
plot(t2, H_all(:,round(position3/dx)))
legendStrings3 = {
    sprintf('n at x = %g um', 10000*position3), ...
    sprintf('m at x = %g um', 10000*position3), ...
    sprintf('h at x = %g um', 10000*position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")
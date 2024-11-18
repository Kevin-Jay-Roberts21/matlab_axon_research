clear all 
close all
clc

dx = 0.0001; % (cm) space step
dt = 0.01; % (ms) time step 
L_my = 0.0075; % (cm) internodal length
L_n = 0.0005; % (cm) nodal length
L_s = L_n + L_my; % (cm) length of an axon segment
n_s = 10; % (dimless) number of axon segments
L = n_s*L_s; % (cm) total length of axon
T = 30; % (ms) the total time of the experiment
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
V_m0 = -64.999; % (mV) initial condition for membrane potential 
V_my0 = -1; % (mV) initial condition for axon potential in periaxonal space
N_0 = 0.3177; % (dimless) initial condition for gating variable n
M_0 = 0.0529; % (dimless) initial condition for gating variable m
H_0 = 0.5961; % (dimless) initial condition for gating variable h

N_n = round(L_n/dx); % number of space steps in a nodal region
N_my = round(L_my/dx); % number of space steps in an internodal region
N_s = N_n + N_my; % number of space steps in an entire axon segement

% defining the total number of space steps and time steps
m = L/dx + 1; % m1 is the number of space steps for V_m
n = T/dt; % n is the number of time steps



% defining the alpha_xi and beta_xi functions
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% defining the b_1(x_i)

B_1 = (a/(2*R_i))*(1 + C_m*a/(C_my*a_my)); % for internodal
B_2 = a/(2*R_i); % for nodal
B_3 = (B_1 + B_2)/2; % for boundary point

% this can be easily verified by testing points like b_1(5) for example
b_1 = @(x) (mod(x - 1, N_s) > N_n).*B_1 + ... % Internodal region
           (mod(x - 1, N_s) < N_n & mod(x - 1, N_s) ~= 0).*B_2 + ... % Nodal region
           ((mod(x - 1, N_s) == N_n) | (mod(x - 1, N_s) == 0)).*((B_1 + B_2)/2); % Boundary point      
      
% defining the initial vectors
Vm = V_m0 * ones(1, m);
Vmy = zeros(1, m);
N = N_0 * ones(1, m);
M = M_0 * ones(1, m);
H = H_0 * ones(1, m);

% putting 0's into the nodal regions of Vmy
for i = 2:m
    if mod(i - 1, N_s) <= N_n
        Vmy(i) = 0; % Nodal region
    else
        Vmy(i) = V_my0; % Internodal region
    end
end

% defining the matrices to collect data at every time step
Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

for j = 1:(n-1)

    %%%%%%%%%%%%%%%%%%%%%%%
    % UPDATING N, M and H %
    %%%%%%%%%%%%%%%%%%%%%%%
    % solving for n^j+1, m^j+1, and h^j+1 given the initial Vm
    for i = 1:m
        newN(i) = 1/(1/dt + alpha_n(Vm(i)) + beta_n(Vm(i))) * (N(i)/dt + alpha_n(Vm(i)));
        newM(i) = 1/(1/dt + alpha_m(Vm(i)) + beta_m(Vm(i))) * (M(i)/dt + alpha_m(Vm(i)));
        newH(i) = 1/(1/dt + alpha_h(Vm(i)) + beta_h(Vm(i))) * (H(i)/dt + alpha_h(Vm(i)));
    end

    % Edit the next U, N, M and H (redefining U, N, M and H vectors)
    N = newN;
    M = newM;
    H = newH;

    % adding all of newN, newM, newH to the _all data
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

    % defining the A matrix and the f vector to solve for V_m
    A = zeros(m, m);
    f = zeros(m, 1);

    % using the boundary conditions to define the top and bottom row of A
    A(1, 1) = 1/dx;
    A(1, 2) = -1/dx;
    A(m, m-1) = -1/dx;
    A(m, m) = 1/dx;

    % creating the first element of V_my, which is 0 since we start in a
    % nodal region
    newVmy(1) = 0;


    % updating the interior of V_m and V_my
    for i = 2:m-1

        % Using the following 3 definitions to determine nodal regions,
        % internodal regions and end points
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment

        % updating the A matrix rows if in an internodal region
        % NOTE: the b function accounts for the piecewise values
        A(i, i-1) = -b_1(i)/dx^2;
        A(i, i) = C_m/dt + 2*b_1(i)/dx^2;
        A(i, i+1) = -b_1(i)/dx^2;

        % If i is at an End Point
        if (i == myelin_start + 1 || i == myelin_end - 1)
            
            % updating the f function at end point: the average (F_1 + F_1)/2
            f(i, 1) = 0.5*((C_m/dt - 1/R_m)*Vm(i) + (1/R_m - C_m/(C_my*R_my))*Vmy(i) + ...
                (C_m/dt - G_K*(newN(i)^4) - G_Na*(newM(i)^3)*newH(i) - G_L)*Vm(i) + ...
                 G_K*(newN(i)^4)*E_L + G_Na*(newM(i)^3)*newH(i)*E_Na + G_L*E_L);
            
            % 0 if in at an end point
            newVmy(i) = 0; 

        % If i is in an Internodal Region
        elseif i >= myelin_start + 1 && i <= myelin_end - 1

            % updating the f function if in an internodal region
            f(i, 1) = (C_m/dt - 1/R_m)*Vm(i) + (1/R_m - C_m/(C_my*R_my))*Vmy(i); 

            % updating Vmy for the internodal region
            newVmy(i) = dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i-1) - ...
                dt*a^2/(C_my*a_my*R_i*dx^2)*Vm(i) + ...
                dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i+1) + ...
                (C_my/dt - 1/R_my)*Vmy(i);   

        % If i is in a Nodal Region
        else 
            
            % updating the f function if in a nodal region
            f(i, 1) = (C_m/dt - G_K*(newN(i)^4) - G_Na*(newM(i)^3)*newH(i) - G_L)*Vm(i) + ...
                G_K*(newN(i)^4)*E_L + G_Na*(newM(i)^3)*newH(i)*E_Na + G_L*E_L; 
            
            % 0 if in a nodal region
            newVmy(i) = 0; 

        end
    end

    j % showing the j index, just for seeing how long simulation takes 

    % because this is technically in a nodal region at the end point
    % recall this experiment is not symmetric
    newVmy(m) = 0;

    % Solving for V_m^{j+1}
    newVm = transpose(A\f);

    % updating Vmy and Vm and adding the data to the _all matrices
    Vm = newVm;
    Vmy = newVmy;

    Vm_all(j+1,:) = Vm;
    Vmy_all(j+1,:) = Vmy;
    N_all(j+1,:) = N;
    M_all(j+1,:) = M;
    H_all(j+1,:) = H;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PICKING TIME AND POSITION SHOTS TO PLOT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

position1 = 0.03; % in cm
position2 = 0.04; % in cm
position3 = 0.05; % in cm
position4 = 0.06; % in cm
position5 = 0.07; % in cm
position6 = 0.075; % in cm
position7 = 0.08; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7];

% Times to observe the voltage along the axon
time1 = 0.5; % in ms
time2 = 1; % in ms
time3 = 1.5; % in ms
time4 = 2; % in ms
time5 = 2.5; % in ms
time6 = 3; % in ms
time7 = 4; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING Vm SPATIAL AND TEMPORAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First figure: Voltage along the axon at different times
figure(1);

% Adjust the figure size (Position [left, bottom, width, height])
set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)

% Create subplot (1 row, 2 columns, 1st subplot)
subplot(1, 2, 1);
t1 = linspace(0, L, m);
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
xlabel("Length of the axon in cm.");

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
    legendStrings2{end+1} = sprintf('$V_m$ at x = %g cm', list_of_positions(i));
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
t1 = linspace(0, L, m);
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
xlabel("Length of the axon in cm.");

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
    legendStrings2{end+1} = sprintf('$V_{my}$ at x = %g cm', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex');
ylabel("$V_{my}$ in millivolts.", 'Interpreter', 'latex');
xlabel("Time in milliseconds.");




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING n, m, and h TEMPROAL PROFILES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(t2, N_all(:,round(position3/dx)))
hold on
plot(t2, M_all(:,round(position3/dx)))
hold on
plot(t2, H_all(:,round(position3/dx)))
legendStrings3 = {
    sprintf('n at x = %g cm', position3), ...
    sprintf('m at x = %g cm', position3), ...
    sprintf('h at x = %g cm', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")
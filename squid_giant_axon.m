clear all
close all
clc

% Tracking the Voltage of a Squid Giant Axon


% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (ohms * cm)
a = 0.025; % axon radius (cm)
L = 5; % axon length (cm)
h = 0.01; % space step (MAY CHANGE LATER)
T = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
g_k = 0.036; % (1/(ohm*cm^2))
g_Na = 0.12; % (1/(ohm*cm^2))
g_L = 0.0003; % (1/(ohm*cm^2))
E_k = -77; % Equilibrium Potential for Potassium Ions (mV)
E_Na = 50; % Equilibrium Potential for Sodium Ions (mV)
E_L = -54.4; % Equilibrium Potential for Leak Channels (mV)

% Mathematical Neuroscience Probability Functions
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% adding sodium conductance (stimulus)
S = 0.2; % (in 1/(ohm*cm^2))
T0 = 5; % start time of when stimulus is added (in ms)
T1 = 5.1; % end time of when stimulus is added (in ms)
P0 = 1; % position of adding the stimulus (in cm)
P1 = 1.1;

% INITIAL CONDITIONS
% for math neuroscience parameters
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
U_0 = -64.9997; % (mV) Voltage (eq: -64.9997)

% number of columns of the matrices (length of axon divided by space step)
m = L/h + 1; 

% initial vectors
U = zeros(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);

% A and b matrix and vector (will be solved in Ax = b style)
A = zeros(m);
b = zeros(m, 1);

% setting the initial U H M and N conditions in the vectors:
U(1,:) = U_0;
N(1,:) = N_0;
M(1,:) = M_0;
H(1,:) = H_0;

% Final matrix with voltage of axon at every time and space step
Uall(1,:) = U;

% Final matrix with probabilities at every time and space step
Nall(1,:) = N;
Mall(1,:) = M; 
Hall(1,:) = H;

% defining new time vectors (each time vector has m elements)
newU = zeros(1, m);
newN = zeros(1, m);
newM = zeros(1, m);
newH = zeros(1, m);

% number of rows in final matrices
n = T/k;

% j is the time step
for j = 1:(n-1)
    
    
    % setting newN, newM, newH vectors (this is a new i, different from the above for loop)
    for i = 1:m
        newN(i) = 1/(1/k + alpha_n(U(i)) + beta_n(U(i))) * (N(i)/k + alpha_n(U(i)));
        newM(i) = 1/(1/k + alpha_m(U(i)) + beta_m(U(i))) * (M(i)/k + alpha_m(U(i)));
        newH(i) = 1/(1/k + alpha_h(U(i)) + beta_h(U(i))) * (H(i)/k + alpha_h(U(i)));
    end
    
    % Edit the next U, N, M and H (redefining U, N, M and H vectors)
    N = newN;
    M = newM;
    H = newH;
    
    % i is the space step
    for i = 1:m

        % defining coefficients
        % Choice 1 set of variables (see powerpoint for Choice 1 & 2 meaning)
        % a1 = -a/(2*r_l*h^2);
        % a2 = a/(r_l*h^2) + c_m/k + g_k*N(i)^4 + g_Na*M(i)^3*H(i) + g_L;
        % a3 = -a/(2*r_l*h^2); 
        % a4 = c_m/k; 
        % a5 = g_k*N(i)^4*E_k + g_Na*M(i)^3*H(i)*E_Na + g_L*E_L;
        
        % Choice 2 set of variables
        a1 = -k*a/(2*r_l*c_m*h^2);
        a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*(N(i)^4)/c_m + k*g_Na*(M(i)^3)*H(i)/c_m + k*g_L/c_m;
        a3 = -k*a/(2*r_l*c_m*h^2); 
        a4 = 1; 
        a5 = k*g_k*N(i)^4*E_k/c_m + k*g_Na*(M(i)^3)*H(i)*E_Na/c_m + k*g_L*E_L/c_m;


        % % adding the stimulus temporally only: (T0 - T1)
        % if j*k >= T0 && j*k <= T1
        %   % must be used for Choice 1
        %   a2 = a/(r_l*h^2) + c_m/k + g_k*N(i)^4 + (g_Na*M(i)^3*H(i) + S) + g_L;
        %   a5 = g_k*N(i)^4*E_k + (g_Na*M(i)^3*H(i) + S)*E_Na + g_L*E_L;
        % 
        %   % must be used for Choice 2
        %   % a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*N(i)^4/c_m + k*(g_Na*M(i)^3*H(i) + S)/c_m + k*g_L/c_m;
        %   % a5 = k*g_k*N(i)^4*E_k/c_m + k*(g_Na*M(i)^3*H(i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
        % end
        
        % adding the stimulus spacially only: (P0 - P1)
        % if i*h >= P0 && i*h <= P1 
        %     % must be used for Choice 1
        %     % a2 = a/(r_l*h^2) + c_m/k + g_k*N(i)^4 + (g_Na*M(i)^3*H(i) + S) + g_L;
        %     % a5 = g_k*N(i)^4*E_k + (g_Na*M(i)^3*H(i) + S)*E_Na + g_L*E_L;
        % 
        %    % must be used for Choice 2
        %    a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*N(i)^4/c_m + k*(g_Na*M(i)^3*H(i) + S)/c_m + k*g_L/c_m;
        %    a5 = k*g_k*N(i)^4*E_k/c_m + k*(g_Na*M(i)^3*H(i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
        % end
        
        
        % adding stimulus temporally and spatially:
        if (j*k >= T0 && j*k <= T1) && (i*h >= P0 && i*h <= P1)
            
            % must be used for Choice 1
            % a2 = a/(r_l*h^2) + c_m/k + g_k*N(i)^4 + (g_Na*M(i)^3*H(i) + S) + g_L;
            % a5 = g_k*N(i)^4*E_k + (g_Na*M(i)^3*H(i) + S)*E_Na + g_L*E_L;
            
            % must be used for Choice 2
            a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*(N(i)^4)/c_m + k*(g_Na*(M(i)^3)*H(i) + S)/c_m + k*g_L/c_m;
            a5 = k*g_k*(N(i)^4)*E_k/c_m + k*(g_Na*(M(i)^3)*H(i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
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

% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = 0.5; % in cm
position2 = 1; % in cm
position3 = 1.5; % in cm
position4 = 2; % in cm
position5 = 2.5; % in cm
position6 = 3; % in cm
position7 = 3.5; % in cm

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
t1 = linspace(0, L, m);
plot(t1, Uall(round(time1/k),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Uall(round(list_of_times(i)/k),:))
end

% describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('Voltage of the axon at time t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex')
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in cm.")

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
    legendStrings2{end+1} = sprintf('Voltage at x = %g cm', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

% plotting N, M, H probability vs time (at a certain position)
figure(3)
plot(t2, Nall(:,round(position3/h)))
hold on
plot(t2, Mall(:,round(position3/h)))
hold on
plot(t2, Hall(:,round(position3/h)))
legendStrings3 = {
    sprintf('N at x = %g cm', position3), ...
    sprintf('M at x = %g cm', position3), ...
    sprintf('H at x = %g cm', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")

% save('salt_cond2023_params_stim.mat');
clear all
close all
clc

% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (ohms * cm)
a = 0.0025; % axon radius (cm)
d = 5; % axon length (cm)
h = 0.01; % space step (MAY CHANGE LATER)
total_time = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
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

% adding sodium conductance (stimulus)
S = 0.003948; % (in 1/(ohm*cm^2))
T0 = 5; % start time of when stimulus is added (in ms)
T1 = 5.1; % end time of when stimulus is added (in ms)
P0 = 1; % position of adding the stimulus (in cm)
P1 = 1.1;


% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

m = d/h; % number of columns of the matrices (length of axon divided by space step)

U = zeros(1, m);
N = zeros(1, m);
M = zeros(1, m);
H = zeros(1, m);
A = zeros(m);
b = zeros(m, 1);
Uall = [];

% setting the initial U H M and N conditions in the vectors:
U(1,:) = V_initial;
N(1,:) = N_0;
M(1,:) = M_0;
H(1,:) = H_0;

Uall(1,:) = U;
Nall(1,:) = N;
Mall(1,:) = M; 
Hall(1,:) = H;

n = total_time/k; % total time is k*n

% j is the time step
for j = 1:(n-1)

    % i is the space step
    for i = 1:m

        % defining coefficients
        a1 = -a/(2*r_l*h^2);
        a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + g_Na*M(1, i)^3*H(1, i) + g_L;
        a3 = -a/(2*r_l*h^2); 
        a4 = c_m/k; 
        a5 = g_k*N(1, i)^4*E_k + g_Na*M(1, i)^3*H(1, i)*E_Na + g_L*E_L;

        % % adding the stimulus during a certain time interval: (T0 - T1)
        % if j*k >= T0 && j*k <= T1
        %     a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
        %     a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % end
        % 
        % adding the stimulus at a spacial interval: (P0 - P1)
        % if i*h >= P0 && i*h <= P1 
        %     a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
        %     a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % end

        % adding stimulus in specific space AND time interval:
        if (j*k >= T0 && j*k <= T1) && (i*h >= P0 && i*h <= P1)
            a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
            a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        end

        % add if statements here for the first row of A and the last row of
        % A
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
            b(i, 1) = a4*U(1, i) + a5; 
        end
    end

    newU = transpose(A\b);
    
    % this is a new i, different from the above for loop
    newN = zeros(1, m);
    newM = zeros(1, m);
    newH = zeros(1, m);
    for i = 1:m
        newN(1, i) = 1/(1/k + alpha_n(U(1, i)) + beta_n(U(1, i))) * (N(1, i)/k + alpha_n(U(1, i)));
        newM(1, i) = 1/(1/k + alpha_m(U(1, i)) + beta_m(U(1, i))) * (M(1, i)/k + alpha_m(U(1, i)));
        newH(1, i) = 1/(1/k + alpha_h(U(1, i)) + beta_h(U(1, i))) * (H(1, i)/k + alpha_h(U(1, i)));
    end

    % Edit the next U, N, M and H (redefining U, N, M and H vectors)
    N(1,:) = newN;
    M(1,:) = newM;
    H(1,:) = newH;
    U(1,:) = newU;

    % USE FOR GRABBING AT EVERY ITERATION
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
% VOLTAGE IS THE SAME AT ANY POSITION
position1 = 0.5; % in cm
position2 = 1; % in cm
position3 = 1.5; % in cm
position4 = 2; % in cm
position5 = 2.5; % in cm
position6 = 3; % in cm
position7 = 3.5; % in cm
position8 = 4; % in cm
position9 = 4.5; % in cm
position10 = 5; % in cm


% Times to observe the voltage along the axon
time1 = 5; % in ms
time2 = 6; % in ms
time3 = 8; % in ms
time4 = 10; % in ms
time5 = 10.5; % in ms
time6 = 10.8; % in ms
time7 = 11; % in ms
time8 = 11.5; % in ms
time9 = 13; % in ms
time10 = 15; % in ms


figure(1)
t1 = linspace(0, d, m);
plot(t1, Uall(time1/k,:))
hold on
plot(t1, Uall(time2/k,:))
hold on
plot(t1, Uall(time3/k,:))
hold on
plot(t1, Uall(time4/k,:))
hold on
plot(t1, Uall(time5/k,:))
hold on
plot(t1, Uall(time6/k,:))
hold on
plot(t1, Uall(time7/k,:))
hold on
plot(t1, Uall(time8/k,:))
hold on
plot(t1, Uall(time9/k,:))
hold on
plot(t1, Uall(time10/k,:))
legend(sprintf('Voltage of the axon at time t = %g ms', time1), sprintf('Voltage of the axon at time t = %g ms', time2), sprintf('Voltage of the axon at time t = %g ms', time3), sprintf('Voltage of the axon at time t = %g ms', time4), sprintf('Voltage of the axon at time t = %g ms', time5), sprintf('Voltage of the axon at time t = %g ms', time6), sprintf('Voltage of the axon at time t = %g ms', time7), sprintf('Voltage of the axon at time t = %g ms', time8), sprintf('Voltage of the axon at time t = %g ms', time9), sprintf('Voltage of the axon at time t = %g ms', time10))
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in cm.")

figure(2)
t2 = linspace(0, total_time, n); % FULL MATRIX
% t2 = linspace(0, total_time, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, total_time, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Uall(:,position1/h))
hold on
plot(t2, Uall(:,position2/h))
hold on
plot(t2, Uall(:,position3/h))
hold on
plot(t2, Uall(:,position4/h))
hold on
plot(t2, Uall(:,position5/h))
hold on
plot(t2, Uall(:,position6/h))
hold on
plot(t2, Uall(:,position7/h))
hold on
plot(t2, Uall(:,position8/h))
hold on
plot(t2, Uall(:,position9/h))
hold on
plot(t2, Uall(:,position10/h))
legend(sprintf('Voltage at x = %g cm', position1),sprintf('Voltage at x = %g cm', position2),sprintf('Voltage at x = %g cm', position3),sprintf('Voltage at x = %g cm', position4),sprintf('Voltage at x = %g cm', position5),sprintf('Voltage at x = %g cm', position6),sprintf('Voltage at x = %g cm', position7),sprintf('Voltage at x = %g cm', position8),sprintf('Voltage at x = %g cm', position9),sprintf('Voltage at x = %g cm', position10))
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

figure(3)
plot(t2, Nall(:,position4/h))
hold on
plot(t2, Mall(:,position4/h))
hold on
plot(t2, Hall(:,position4/h))
legend(sprintf('N at x = %g cm', position4), sprintf('M at x = %g cm', position4), sprintf('H at x = %g cm', position4))
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")





% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )

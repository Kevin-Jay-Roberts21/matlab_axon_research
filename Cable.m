clear all
close all
clc

% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (mS / (ohm*cm^2))
r_l = 30; % membrane resistance (ohms * cm)
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

% adding sodium conductance (sitmulus)
S = 0.01;
T0 = 5; % start time of when stimulus is added (in ms)
T1 = 5; % end time of when stimulus is added (in ms)

% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

m = d/h; % number of columns of the matrices (length of axon divided by space step)

U = zeros(m, 1);
N = zeros(m, 1);
M = zeros(m, 1);
H = zeros(m, 1);
A = zeros(m);
b = zeros(m, 1);
Uall = [];

% setting the initial U H M and N conditions in the vectors:
U(:,1) = V_initial;
N(:,1) = N_0;
M(:,1) = M_0;
H(:,1) = H_0;

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
        a2 = a/(r_l*h^2) + c_m/k + g_k*N(i, 1)^4 + g_Na*M(i, 1)^3*H(i, 1) + g_L;
        a3 = -a/(2*r_l*h^2); 
        a4 = c_m/k; 
        a5 = g_k*N(i, 1)^4*E_k + g_Na*M(i, 1)^3*H(i, 1)*E_Na + g_L*E_L;

        % adding the stimulus at 5 ms.
        if j*k >= T0 && j*k <= T1
            a2 = a/(r_l*h^2) + c_m/k + g_k*N(i, 1)^4 + (g_Na*M(i, 1)^3*H(i, 1) + S) + g_L;
            a5 = g_k*N(i, 1)^4*E_k + (g_Na*M(i, 1)^3*H(i, 1) + S)*E_Na + g_L*E_L;
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
            b(i, 1) = a4*U(i, 1) + a5; 
        end
    end

    newU = A\b;

    newN = 1/(1/k + 0.01*(U(i, 1) + 55)/(1 - exp(-(U(i, 1) + 55)/10)) + 0.125*exp(-(U(i, 1) + 65)/80)) * (N(i, 1)/k + 0.01*(U(i, 1) + 55)/(1 - exp(-(U(i, 1) + 55)/10)));
    newM = 1/(1/k + 0.1*(U(i, 1) + 40)/(1 - exp(-(U(i, 1) + 40)/10)) + 4*exp(-(U(i, 1) + 65)/18)) * (M(i, 1)/k + 0.1*(U(i, 1) + 40)/(1 - exp(-(U(i, 1) + 40)/10)));
    newH = 1/(1/k + 0.07*exp(-(U(i, 1) + 65)/20) + 1/(1 + exp(-(U(i, 1) + 35)/10))) * (H(i, 1)/k + 0.07*exp(-(U(i, 1) + 65)/20));

    % Edit the next U, N, M and H (redefining U, N, M and H vectors)
    N(:,1) = newN; 
    M(:,1) = newM;
    H(:,1) = newH;
    U(:,1) = newU;

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

figure(1)
t1 = linspace(0, 5, m);
plot(t1, U)
legend("Voltage at end of time across the axon.")

figure(2)
% now pick a position to plot all of the voltages
% VOLTAGE IS THE SAME AT ANY POSITION
position1 = 30; % (this is at positon x = 30*h cm = 30*0.01 = 0.3cm)
t2 = linspace(0, total_time, n); % FULL MATRIX
% t2 = linspace(0, total_time, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, total_time, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Uall(:,position1))
legend("Voltage at x = 0.3cm")
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

figure(3)
plot(t2, Nall(:,position1))
hold on
plot(t2, Mall(:,position1))
hold on
plot(t2, Hall(:,position1))
legend("N at x = 0.3cm", "M at x = 0.3cm", "H at x = 0.3cm")
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")





% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )
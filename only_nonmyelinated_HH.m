clear all
close all
clc

% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (ohms * cm)
a = 0.0025; % axon radius (cm)
d = 1.0005; % axon length (cm)
h = 0.0005; % space step (MAY CHANGE LATER)
T = 20; % we only ever want to run up to 35 ms (where we find equilibrium)
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
S = 0.234148; % (in 1/(ohm*cm^2))
T0 = 2; % start time of when stimulus is added (in ms)
T1 = 2.1; % end time of when stimulus is added (in ms)
P0 = 0.1000; % position of adding the stimulus (in cm)
P1 = 0.1005;


% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

m = round(d/h); % number of columns of the matrices (length of axon divided by space step)

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

n = T/k; % total time is k*n

tic; % tracking the time

% j is the time step
for j = 1:(n-1)

    % i is the space step
    for i = 1:m

        % defining coefficients
        % first set of a's
        % a1 = -a/(2*r_l*h^2);
        % a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + g_Na*M(1, i)^3*H(1, i) + g_L;
        % a3 = -a/(2*r_l*h^2); 
        % a4 = c_m/k; 
        % a5 = g_k*N(1, i)^4*E_k + g_Na*M(1, i)^3*H(1, i)*E_Na + g_L*E_L;
        % second set of a's
        a1 = -k*a/(2*r_l*c_m*h^2);
        a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*(N(1, i)^4)/c_m + k*g_Na*(M(1, i)^3)*H(1, i)/c_m + k*g_L/c_m;
        a3 = -k*a/(2*r_l*c_m*h^2); 
        a4 = 1; 
        a5 = k*g_k*N(1, i)^4*E_k/c_m + k*g_Na*(M(1, i)^3)*H(1, i)*E_Na/c_m + k*g_L*E_L/c_m;



        % % adding the stimulus during a certain time interval: (T0 - T1)
        % if j*k >= T0 && j*k <= T1
        %   % must be used for the first set of a's
        %   a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
        %   a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % 
        %   % must be used for the second set of a's
        %   % a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*N(1, i)^4/c_m + k*(g_Na*M(1, i)^3*H(1, i) + S)/c_m + k*g_L/c_m;
        %   % a5 = k*g_k*N(1, i)^4*E_k/c_m + k*(g_Na*M(1, i)^3*H(1, i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
        % end
        % 
        % adding the stimulus at a spacial interval: (P0 - P1)
        % if i*h >= P0 && i*h <= P1 
        %     % must be used for the first set of a's
        %     % a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
        %     % a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % 
        %    % must be used for the second set of a's
        %    a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*N(1, i)^4/c_m + k*(g_Na*M(1, i)^3*H(1, i) + S)/c_m + k*g_L/c_m;
        %    a5 = k*g_k*N(1, i)^4*E_k/c_m + k*(g_Na*M(1, i)^3*H(1, i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
        % end
        % adding stimulus in specific space AND time interval:
        if (j*k >= T0 && j*k <= T1) && (i*h >= P0 && i*h <= P1)
            
            % must be used for the first set of a's
            % a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
            % a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
            
            % must be used for the second set of a's
            a2 = 1 + k*a/(r_l*c_m*h^2) + k*g_k*(N(1, i)^4)/c_m + k*(g_Na*(M(1, i)^3)*H(1, i) + S)/c_m + k*g_L/c_m;
            a5 = k*g_k*(N(1, i)^4)*E_k/c_m + k*(g_Na*(M(1, i)^3)*H(1, i) + S)*E_Na/c_m + k*g_L*E_L/c_m;
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

elapsedtime = toc % getting the end time
max(Uall(:))

% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = 0.001; % in cm 
position2 = 0.0015;
position3 = 0.02; 
position4 = 0.03; 
position5 = 0.1; 
position6 = 0.2; 
position7 = 0.3; 
position8 = 0.4;
position9 = 0.5;
position10 = 0.6;
% position1 = 0.5; % in cm
% position2 = 1; % in cm
% position3 = 1.5; % in cm
% position4 = 2; % in cm
% position5 = 2.5; % in cm
% position6 = 3; % in cm
% position7 = 3.5; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7];

% Times to observe the voltage along the axon
time1 = 2; % in ms
time2 = 2.1; % in ms
time3 = 5; % in ms
time4 = 10; % in ms
time5 = 11; % in ms
time6 = 11.5; % in ms
time7 = 12; % in ms
% time1 = 5; % in ms
% time2 = 6; % in ms
% time3 = 8; % in ms
% time4 = 10; % in ms
% time5 = 10.5; % in ms
% time6 = 10.8; % in ms
% time7 = 11; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];

figure(1)
t1 = linspace(0, d, m);
legendStrings1 = {};
plot(t1, Uall(round(time1/k),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Uall(round(list_of_times(i)/k),:))
end
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('Voltage of the axon at time t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex')
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in cm.")

figure(2)
t2 = linspace(0, T, n); % FULL MATRIX
% t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
legendStrings2 = {};
plot(t2, Uall(:,round(position1/h)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, Uall(:,round(list_of_positions(i)/h)))
end
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('Voltage at x = %g cm', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

% printing the repolarization information
% for i = 1:length(list_of_positions)
%     [speed, time_difference, voltage_difference] = repolarization_function(Uall, list_of_positions(i), V_initial, h, k)
% end


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





% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )

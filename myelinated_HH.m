clear all
close all
clc

% SPECIAL NOTE: to get results in terms of micrometers, just change the 
% length of the axon to 0.1 for example and shrink the space step.

% defining all of the initial and constant variables
c_m_nodal = 0.002; % membrane capacitance (ms / (ohm*cm^2))
c_m_internodal = 0.001; % (ms / (ohm*cm^2))
r_l = 70; % specific intracellular resistivity (or axoplasmic resistivity) (ohms * cm)
radius_nodal = 0.000165; % (cm)
radius_internodal = 0.0003; % (cm)
h = 0.0001; % space step (this)
total_time = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
g_L = 0.08; % (1/(ohm*cm^2))
g_k_nodal = 0.08; % (1/(ohm*cm^2))
g_k_internodal = 0; % (1/(ohm*cm^2))
g_Na_nodal = 3; % (1/(ohm*cm^2))
g_Na_internodal = 0; % (1/(ohm*cm^2))
E_k = -84; % (mV)
E_Na = 50; % (mV)
E_L = -83.38; % (mV)

% The following are functions of voltage
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% defining nodal regions, and the axon length will be based on how many
% regions we have 
num_of_nodes = 7;
% the first number is in um, the *10^(-4) converts it to cm
nodal_length = 5*10^(-4); % (in cm)
myelinated_length = 95*10^(-4); % (in cm)
d = (nodal_length * num_of_nodes) + (myelinated_length * (num_of_nodes - 1)); % axon length (in cm)


% creating a list of nodel regions [[start_pos1, end_pos1], [start_pos2, end_pos2], ...]
nodal_regions = [];
for i = 0:(num_of_nodes-1)
    nodal_regions(:,i+1) = [(i*nodal_length)+(i*myelinated_length), (i*nodal_length)+(i*myelinated_length)+nodal_length];
end

% (the following conductances and radius are fcns of space, to mimmick myelin)
% starting each of the functions with a condition. This condition being the
% values of each if the position is in the first nodal region
a = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
c_m = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
g_Na = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));
g_k = @(x) (x >= nodal_regions(1, 1) && x <= nodal_regions(2, 1));

% creating a for loop that will turn the a, g_k and g_Na into piecewise
% functions, having different values in nodal regions vs myelinated regions 
for i = 2:num_of_nodes
    a = @(x) a(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    c_m = @(x) c_m(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    g_Na = @(x) g_Na(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
    g_k = @(x) g_k(x) || (x >= nodal_regions(1, i) && x <= nodal_regions(2, i));
end

% finally, we make a, g_Na, and g_k certain values if the conditions above
% are met or not the first element of the summation is if x is in a nodal
% region, the second element is if the x is in a myelinated region.
a = @(x) radius_nodal*a(x) + radius_internodal*(a(x)==0); % axon radius (where both values (0.0025 and 0.005) are in cm)
c_m = @(x) c_m_nodal*c_m(x) + c_m_internodal*(c_m(x)==0);
g_k = @(x) g_k_nodal*g_k(x) + g_k_internodal*(g_k(x)==0); % (g_k is in 1/(ohm*cm^2))
g_Na = @(x) g_Na_nodal*g_Na(x) + g_Na_internodal*(g_Na(x)==0); % (g_Na is in 1/(ohm*cm^2))


% adding sodium conductance (stimulus)
S = 10; %0.01; % (in 1/(ohm*cm^2))
T0 = 5; % start time of when stimulus is added (in ms)
T1 = 6; % end time of when stimulus is added (in ms)

% NOTE: the stimulus MUST be added in a nodal region from (0um to 1um is fine)
P0 = 0; % position of adding the stimulus (in cm)
P1 = 0.0002; % ending position of adding the stimulus (in cm or *10^4 in um)

% INITIAL CONDITIONS
N_0 = 0.1009; % probability that potassium gate is open (eq: 0.1009)
M_0 = 0.0051; % probability that Sodium activation gate is open (eq: 0.0051)
H_0 = 0.9571; % probability that Sodium inactivation gate is open (eq: 0.9571)
V_initial = -83.3794; % (mV) Voltage (eq: -83.3794)

% should be CAREFUL about rounding here. For some reason matlab thinks that
% something like 3.678e3 is not considered an interger sometimes
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

n = total_time/k; % total time is k*n

% j is the time step
for j = 1:(n-1)

    % i is the space step
    for i = 1:m

        % defining coefficients
        a1 = -a(i*h)/(2*r_l*h^2);
        a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + g_Na(i*h)*M(1, i)^3*H(1, i) + g_L;
        a3 = -a(i*h)/(2*r_l*h^2); 
        a4 = c_m(i*h)/k; 
        a5 = g_k(i*h)*N(1, i)^4*E_k + g_Na(i*h)*M(1, i)^3*H(1, i)*E_Na + g_L*E_L;

        % % adding the stimulus during a certain time interval: (T0 - T1)
        % if j*k >= T0 && j*k <= T1
        %     a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
        %     a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % end
        % 
        % % adding the stimulus at a spacial interval: (P0 - P1)
        % if i*h >= P0 && i*h <= P1 
        %     a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
        %     a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % end

        % % adding stimulus in specific space AND time interval:
        if (j*k >= T0 && j*k <= T1) && (i*h >= P0 && i*h <= P1)
            a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
            a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
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

% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = 0.0001; % in cm 
position2 = 0.0050; % in cm
position3 = 0.0100; % in cm
position4 = 0.0150; % in cm
position5 = 0.0200; % in cm 
position6 = 0.0250; % in cm
position7 = 0.0350; % in cm
position8 = 0.0450; % in cm
position9 = 0.0550; % in cm 
position10 = 0.0600; % in cm

% Times to observe the voltage along the axon
time1 = 1; % in ms
time2 = 5; % in ms
time3 = 5.5; % in ms
time4 = 6; % in ms
time5 = 6.2; % in ms
time6 = 6.5; % in ms
time7 = 7; % in ms
time8 = 8; % in ms
time9 = 9; % in ms
time10 = 10; % in ms


d_in_um = round(d*10000); % using this value to display plots in um instead of cm

figure(1)
t1 = linspace(0, d_in_um, m);
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
xlabel("Length of the axon in um.")

figure(2)
t2 = linspace(0, total_time, n); % FULL MATRIX
% t2 = linspace(0, total_time, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, total_time, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Uall(:,round(position1/h)))
hold on
plot(t2, Uall(:,round(position2/h)))
hold on
plot(t2, Uall(:,round(position3/h)))
hold on
plot(t2, Uall(:,round(position4/h)))
hold on
plot(t2, Uall(:,round(position5/h)))
hold on
plot(t2, Uall(:,round(position6/h)))
hold on
plot(t2, Uall(:,round(position7/h)))
hold on
plot(t2, Uall(:,round(position8/h)))
hold on
plot(t2, Uall(:,round(position9/h)))
hold on
plot(t2, Uall(:,round(position10/h)))
legend(sprintf('Voltage at x = %g um', position1*10000),sprintf('Voltage at x = %g um', position2*10000),sprintf('Voltage at x = %g um', position3*10000),sprintf('Voltage at x = %g um', position4*10000),sprintf('Voltage at x = %g um', position5*10000),sprintf('Voltage at x = %g um', position6*10000),sprintf('Voltage at x = %g um', position7*10000),sprintf('Voltage at x = %g um', position8*10000),sprintf('Voltage at x = %g um', position9*10000),sprintf('Voltage at x = %g um', position10*10000))
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

figure(3)
plot(t2, Nall(:,round(position4/h)))
hold on
plot(t2, Mall(:,round(position4/h)))
hold on
plot(t2, Hall(:,round(position4/h)))
legend(sprintf('N at x = %g cm', position4), sprintf('M at x = %g cm', position4), sprintf('H at x = %g cm', position4))
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")


% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )

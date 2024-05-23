clear all
close all
clc

% in this file, I use a radius and length of axon in terms of micrometers
% instead of centimeters. We don't have to convert anything else to
% micrometers as these lengths will cancel themselves out.

% In this file, as opposed to having everything in terms of cm, I converted
% everything to be in terms of micrometers or um. NOTE: WE MUST take into
% account that this will change the space step, we may now make the space
% step quite a bit larger than h = 0.01.

% defining all of the initial and constant variables
c_m = 0.001; % membrane capacitance (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (ohms * cm)
a = 25; % axon radius (um)
h = 0.1; % space step (MAY CHANGE LATER)
total_time = 35; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
g_Na = 0.12; % (1/(ohm*cm^2))
g_k = 0.036; % (1/(ohm*cm^2))
g_L = 0.0003; % (1/(ohm*cm^2))
E_Na = 50; % (mV)
E_k = -77; % (mV)
E_L = -54.4; % (mV)

% The following are functions of voltage
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% adding sodium conductance (stimulus)
S = 0.00017; % (in 1/(ohm*um^2))
T0 = 0; % start time of when stimulus is added (in ms)
T1 = 0.1; % end time of when stimulus is added (in ms)

% NOTE: the stimulus MUST be added in a nodal region from (0um to 1um is fine)
P0 = 0; % position of adding the stimulus (in um)
P1 = 1; % ending position of adding the stimulus (in um)

% defining nodal regions, and the axon length will be based on how many
% regions we have
num_of_nodes = 3
nodal_length = 2.3 % (in um)
myelinated_length = 85 % (in um)
d = (nodal_length * num_of_nodes) + (myelinated_length * num_of_nodes) % axon length (in cm)

% creating a list of nodel regions [[start_pos1, end_pos1], [start_pos2, end_pos2], ...]

nodal_regions = []
for i = 0:num_of_nodes
    nodal_regions(:,i+1) = [(i*nodal_length)+(i*myelinated_length), (i*nodal_length)+(i*myelinated_length)+nodal_length];
end

nodal_regions(2, 1)


% INITIAL CONDITIONS
N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

% should be careful about rounding here. For some reason matlab thinks that
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
        
        % using l as the index   
        for l = 1:size(nodal_regions, 2) % for the number of columns in nodal_regions

            % if inside a nodal region
            if i*h >= nodal_regions(1, l) && i*h <= nodal_regions(2, l)
                g_k = g_k; 
                g_Na = g_Na;
                g_L = g_L;

            % if inside a myelinated region
            else
                g_k = 0; 
                g_Na = 0;
                g_L = 0;
            end
        end 

        % defining coefficients
        a1 = -a/(2*r_l*h^2);
        a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + g_Na*M(1, i)^3*H(1, i) + g_L;
        a3 = -a/(2*r_l*h^2); 
        a4 = c_m/k; 
        a5 = g_k*N(1, i)^4*E_k + g_Na*M(1, i)^3*H(1, i)*E_Na + g_L*E_L;

        % adding the stimulus during a certain time interval: (T0 - T1)
        % if i*h >= T0 && i*h <= T1
        %     a2 = a/(r_l*h^2) + c_m/k + g_k*N(1, i)^4 + (g_Na*M(1, i)^3*H(1, i) + S) + g_L;
        %     a5 = g_k*N(1, i)^4*E_k + (g_Na*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % end
            
        


        % adding the stimulus every in the second nodal section
        if i*h >= nodal_regions(1, 2) && i*h <= nodal_regions(2, 2)
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
position1 = 1; % in um
position2 = 100; % in um
position3 = 200; % in um
position4 = 250; % in um



% Times to observe the voltage along the axon
time1 = 5; % in ms
time2 = 8; % in ms
time3 = 12; % in ms
time4 = 14; % in ms

figure(1)
t1 = linspace(0, d, m);
plot(t1, Uall(time1/h,:))
hold on
plot(t1, Uall(time2/h,:))
hold on
plot(t1, Uall(time3/h,:))
hold on
plot(t1, Uall(time4/h,:))
legend(sprintf('Voltage of the axon at time t = %g ms', time1), sprintf('Voltage of the axon at time t = %g ms', time2), sprintf('Voltage of the axon at time t = %g ms', time3), sprintf('Voltage of the axon at time t = %g ms', time4))
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in um.")

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
legend(sprintf('Voltage at x = %g um', position1),sprintf('Voltage at x = %g um', position2),sprintf('Voltage at x = %g um', position3),sprintf('Voltage at x = %g um', position4))
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

figure(3)
plot(t2, Nall(:,position1/h))
hold on
plot(t2, Mall(:,position1/h))
hold on
plot(t2, Hall(:,position1/h))
legend(sprintf('N at x = %g um', position4), sprintf('M at x = %g um', position4), sprintf('H at x = %g um', position4))
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")





% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )

clear all
close all
clc

% SPECIAL NOTES: to get results in terms of micrometers, just change the 
% length of the axon to 0.1 for example and shrink the space step.

% defining all of the initial and constant variables
c_m_nodal = 0.001; % membrane capacitance (ms / (ohm*cm^2))
c_m_internodal = 0.001; % (ms / (ohm*cm^2))
r_l = 30; % specific intracellular resistivity (or axoplasmic resistivity) (ohms * cm)
radius_nodal = 0.0025; % (cm)
radius_internodal = 0.0025; % (cm)
h = 0.0005; % space step (this)
T = 15; % we only ever want to run up to 35 ms (where we find equilibrium)
k = 0.01; % time step (MAY CHANGE LATER)
g_L = 0.0003; % (1/(ohm*cm^2))
g_k_nodal = 0.036; % (1/(ohm*cm^2))
g_k_internodal = 0; % (1/(ohm*cm^2))
g_Na_nodal = 0.12; % (1/(ohm*cm^2))
g_Na_internodal = 0; % (1/(ohm*cm^2))
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

% defining nodal regions, and the axon length will be based on how many
% regions we have 
num_of_nodes = 31;

% make the ion channels uniform



% the first number is in um, the *10^(-4) converts it to cm
nodal_length = 0.0010; % (in cm)
myelinated_length = 0.0115; % (in cm)
d = (myelinated_length*(num_of_nodes-1)) + (nodal_length*(num_of_nodes-1)) + nodal_length; % axon length (in cm)


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
S = 10; % (in 1/(ohm*cm^2))
T0 = 2; % start time of when stimulus is added (in ms)
T1 = 2.1; % end time of when stimulus is added (in ms)
% NOTE: the stimulus MUST be added in a nodal region
P0 = 0.1000; % position of adding the stimulus (in cm)
P1 = 0.1005; % ending position of adding the stimulus (in cm or *10^4 in um)

% INITIAL CONDITIONS
% N_0 = 0.4749; % probability that potassium gate is open (eq: 0.3177)
% M_0 = 0.1575; % probability that Sodium activation gate is open (eq: 0.0529)
% H_0 = 0.2636; % probability that Sodium inactivation gate is open (eq: 0.5961)
% V_initial = -55.4388; % (mV) Voltage (eq: -64.9997)
N_0 = 0.4672; % probability that potassium gate is open (eq: 0.3177)
M_0 = 0.1499; % probability that Sodium activation gate is open (eq: 0.0529)
H_0 = 0.2770; % probability that Sodium inactivation gate is open (eq: 0.5961)
V_initial = -55.5340; % (mV) Voltage (eq: -64.9997)
% N_0 = 0.3177; % probability that potassium gate is open (eq: 0.3177)
% M_0 = 0.0529; % probability that Sodium activation gate is open (eq: 0.0529)
% H_0 = 0.5961; % probability that Sodium inactivation gate is open (eq: 0.5961)
% V_initial = -64.9997; % (mV) Voltage (eq: -64.9997)

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

n = T/k; % total time is k*n

% j is the time step
for j = 1:(n-1)
    j
    % i is the space step
    for i = 1:m

        % defining coefficients
        % first set of a's
        % a1 = -a(i*h)/(2*r_l*h^2);
        % a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + g_Na(i*h)*M(1, i)^3*H(1, i) + g_L;
        % a3 = -a(i*h)/(2*r_l*h^2); 
        % a4 = c_m(i*h)/k; 
        % a5 = g_k(i*h)*N(1, i)^4*E_k + g_Na(i*h)*M(1, i)^3*H(1, i)*E_Na + g_L*E_L;
        % second set of a's
        a1 = -k*a(i*h)/(2*r_l*c_m(i*h)*h^2);
        a2 = 1 + k*a(i*h)/(r_l*c_m(i*h)*h^2) + k*g_k(i*h)*(N(1, i)^4)/c_m(i*h) + k*g_Na(i*h)*(M(1, i)^3)*H(1, i)/c_m(i*h) + k*g_L/c_m(i*h);
        a3 = -k*a(i*h)/(2*r_l*c_m(i*h)*h^2); 
        a4 = 1; 
        a5 = k*g_k(i*h)*(N(1, i)^4)*E_k/c_m(i*h) + k*g_Na(i*h)*(M(1, i)^3)*H(1, i)*E_Na/c_m(i*h) + k*g_L*E_L/c_m(i*h);



        % % adding the stimulus during a certain time interval: (T0 - T1)
        % if j*k >= T0 && j*k <= T1
        %   % must be used for the first set of a's
        %   a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
        %   a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % 
        %   % must be used for the second set of a's
        %   % a2 = 1 + k*a(i*h)/(r_l*c_m(i*h)*h^2) + k*g_k(i*h)*N(1, i)^4/c_m(i*h) + k*(g_Na(i*h)*M(1, i)^3*H(1, i) + S)/c_m(i*h) + k*g_L/c_m(i*h);
        %   % a5 = k*g_k(i*h)*N(1, i)^4*E_k/c_m(i*h) + k*(g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na/c_m(i*h) + k*g_L*E_L/c_m(i*h);
        % end
        % 
        % adding the stimulus at a spacial interval: (P0 - P1)
        % if i*h >= P0 && i*h <= P1 
        %     % must be used for the first set of a's
        %     % a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
        %     % a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
        % 
        %    % must be used for the second set of a's
        %    a2 = 1 + k*a(i*h)/(r_l*c_m(i*h)*h^2) + k*g_k(i*h)*N(1, i)^4/c_m(i*h) + k*(g_Na(i*h)*M(1, i)^3*H(1, i) + S)/c_m(i*h) + k*g_L/c_m(i*h);
        %    a5 = k*g_k(i*h)*N(1, i)^4*E_k/c_m(i*h) + k*(g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na/c_m(i*h) + k*g_L*E_L/c_m(i*h);
        % end
        % adding stimulus in specific space AND time interval:
        if (j*k >= T0 && j*k <= T1) && (i*h >= P0 && i*h <= P1)
            
            % must be used for the first set of a's
            % a2 = a(i*h)/(r_l*h^2) + c_m(i*h)/k + g_k(i*h)*N(1, i)^4 + (g_Na(i*h)*M(1, i)^3*H(1, i) + S) + g_L;
            % a5 = g_k(i*h)*N(1, i)^4*E_k + (g_Na(i*h)*M(1, i)^3*H(1, i) + S)*E_Na + g_L*E_L;
            
            % must be used for the second set of a's
            a2 = 1 + k*a(i*h)/(r_l*c_m(i*h)*h^2) + k*g_k(i*h)*(N(1, i)^4)/c_m(i*h) + k*(g_Na(i*h)*(M(1, i)^3)*H(1, i) + S)/c_m(i*h) + k*g_L/c_m(i*h);
            a5 = k*g_k(i*h)*(N(1, i)^4)*E_k/c_m(i*h) + k*(g_Na(i*h)*(M(1, i)^3)*H(1, i) + S)*E_Na/c_m(i*h) + k*g_L*E_L/c_m(i*h);
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
time3 = 3; % in ms
time4 = 4; % in ms
time5 = 5; % in ms
time6 = 6; % in ms
time7 = 7; % in ms

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
plot(t2, Nall(:,round(position1/h)))
hold on
plot(t2, Mall(:,round(position1/h)))
hold on
plot(t2, Hall(:,round(position1/h)))
legendStrings3 = {
    sprintf('N at x = %g cm', position3), ...
    sprintf('M at x = %g cm', position3), ...
    sprintf('H at x = %g cm', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")


% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )

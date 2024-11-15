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

N_n = round(L_n / dx); % number of space steps in a nodal region
N_my = round(L_my / dx); % number of space steps in an internodal region
N_s = N_n + N_my; % number of space steps in an entire axon segement

% defining the total number of space steps and time steps
m1 = L/dx + 1; % m1 is the number of space steps for V_m
m2 = N_my * n_s; % m2 is the number of space steps for V_my
n = T/dt; % n is the number of time steps



% defining the alpha_xi and beta_xi functions
alpha_n = @(V) 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
beta_n = @(V) 0.125*exp(-(V + 65)/80);
alpha_m = @(V) 0.1*(V + 40)/(1 - exp(-(V + 40)/10));
beta_m = @(V) 4*exp(-(V + 65)/18);
alpha_h = @(V) 0.07*exp(-(V + 65)/20);
beta_h = @(V) 1/(1 + exp(-(V + 35)/10));

% defining the b_1(x_i), f_1(x_i), and f_2(x_i)
b_1 = @(x) (mod(x, N_s) > N_n).*(a/(2*R_i)).*(1 + C_m*a/(C_my*a_my)) + ...
          (mod(x, N_s) <= N_n).*(a/(2*R_i));

% defining the initial vectors
Vm = V_m0 * ones(1, m1);
Vmy = V_my0 * ones(1, m1);
N = N_0 * ones(1, m1);
M = M_0 * ones(1, m1);
H = H_0 * ones(1, m1);

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
    for i = 1:m1
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
    A = zeros(m1, m1);
    f = zeros(m1, 1);
    
    % using the boundary conditions to define the top and bottom row of A
    A(1, 1) = 1/dx;
    A(1, 2) = -1/dx;
    A(m1, m1-1) = -1/dx;
    A(m1, m1) = 1/dx;
    
    % creating the first element of V_my, which is 0 since we start in a
    % nodal region
    newVmy(1) = 0;
    

    % updating the interior of V_m and V_my
    for i = 2:m1-1

        
        % solving for Vmy^(j+1) using Vm^j
        % Check if index i is in an internodal region
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % Start of internodal region in this segment
        myelin_end = seg*(N_s); % End of internodal region in this segment
        
       
        % updating the A matrix rows if in an internodal region
        % NOTE: the b function accounts for the piecewise values
        A(i, i-1) = -b_1(i)/dx^2;
        A(i, i) = C_m/dt + 2*b_1(i)/dx^2;
        A(i, i+1) = -b_1(i)/dx^2;
        
        
        % Condition to check if i is within an internodal region
        if i >= myelin_start + 1 && i <= myelin_end - 1
            
            % updating the f function if in an internodal region
            f(i, 1) = (C_m/dt - 1/R_m)*Vm(i) + (1/R_m - C_m/(C_my*R_my))*Vmy(i); 
            
            %%%%%%%%%%%%%%%%
            % UPDATING Vmy %
            %%%%%%%%%%%%%%%%
            newVmy(i) = dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i-1) - ...
                dt*a^2/(C_my*a_my*R_i*dx^2)*Vm(i) + ...
                dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(i+1) + ...
                (C_my/dt - 1/R_my)*Vmy(i);   
        
        % within a nodal region
        else 
            
            newVmy(i) = 0; % 0 if in a nodal region
            
            % updating the f function if in a nodal region
            f(i, 1) = (C_m/dt - G_K*(newN(i)^4) - G_Na*(newM(i)^3)*newH(i) - G_L)*Vm(i) + ...
                G_K*(newN(i)^4)*E_L + G_Na*(newM(i)^3)*newH(i)*E_Na + G_L*E_L; 
        end
    end
    
    % becuase this is technically in a nodal region
    newVmy(m1) = 0;

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

% plotting Voltage vs Axon length
figure(1)
t1 = linspace(0, L, m1);
plot(t1, Vm_all(round(time1/dt),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Vm_all(round(list_of_times(i)/dt),:))
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
length(t2)
length(Vm_all(:,round(position1/dx)))

plot(t2, Vm_all(:,round(position1/dx)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, Vm_all(:,round(list_of_positions(i)/dx)))
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
plot(t2, N_all(:,round(position3/dx)))
hold on
plot(t2, M_all(:,round(position3/dx)))
hold on
plot(t2, H_all(:,round(position3/dx)))
legendStrings3 = {
    sprintf('N at x = %g cm', position3), ...
    sprintf('M at x = %g cm', position3), ...
    sprintf('H at x = %g cm', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")
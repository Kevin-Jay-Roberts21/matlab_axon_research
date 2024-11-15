clear all 
close all
clc

dx = 0.0001; % (cm) space step
dt = 0.01; % (ms) time step 
L_my = 0.0075; % (cm) internodal length
L_n = 0.0005; % (cm) nodal length
L_s = L_n + L_my; % (cm) length of an axon segment
n_s = 5; % (dimless) number of axon segments
L = n_s*L_s; % (cm) total length of axon
T = 5; % (ms) the total time of the experiment
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
Vmy = V_my0 * ones(1, m2);
N = N_0 * ones(1, m1);
M = M_0 * ones(1, m1);
H = H_0 * ones(1, m1);

Vm_all(1,:) = Vm;
Vmy_all(1,:) = Vmy;
N_all(1,:) = N;
M_all(1,:) = M; 
H_all(1,:) = H;

length(Vm)
length(Vmy)

for j = 1:n

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
            % updating the f function if in a nodal region
            f(i, 1) = (C_m/dt - G_K*(newN(i)^4) - G_Na*(newM(i)^3)*newH(i) - G_L)*Vm(i) + ...
                G_K*(newN(i)^4)*E_L + G_Na*(newM(i)^3)*newH(i)*E_Na + G_L*E_L; 
        end
    end
    
    % handling the neumann boundary conditions for Vm 
    % (i.e. (V_{m,m1} - V_{m,m1-1})/dx = 0 => V_{m,m1} = V_{m,m1-1})
    newVmy(m2) = dt*a^2/(2*C_my*a_my*R_i*dx^2)*Vm(m2-1) - ...
        2*dt*a^2/(C_my*a_my*R_i*dx^2)*Vm(m2) + ...
        (1 - dt/(C_my*R_my))*Vmy(m2);  

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

length(newVm)
length(newVmy)
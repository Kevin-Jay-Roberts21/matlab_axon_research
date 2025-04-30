% HH model
% Kevin Roberts
% April 2025

% This function will solve for Vm, n, m and h for the HH model given mesh
% and material parameters

function hh_solution_data = hh_function_v1(mesh_params, material_params)
    
    % Grabing and defining all the inputed mesh and material parameters
    
    % mesh params
    L = mesh_params.L;
    dx = mesh_params.dx;
    T = mesh_params.T;
    dt = mesh_params.dt;
    m = mesh_params.m;
    n = mesh_params.n;
    
    % material params
    a = material_params.a;
    C_m = material_params.C_m;
    R_i = material_params.R_i;
    G_K = material_params.G_K;
    G_Na = material_params.G_Na;
    G_L = material_params.G_L;
    E_K = material_params.E_K;
    E_Na = material_params.E_Na;
    E_L = material_params.E_L;
    T_base = material_params.T_base;
    T_actual = material_params.T_actual;
    Q_10 = material_params.Q_10;
    S_v = material_params.S_v;
    S_T0 = material_params.S_T0; 
    S_T1 = material_params.S_T1; 
    S_P0 = material_params.S_P0; 
    S_P1 = material_params.S_P1; 
    N_0 = material_params.N_0;
    M_0 = material_params.M_0; 
    H_0 = material_params.H_0;
    V_m0 = material_params.V_m0;

    % Defining the temperature factor and the alpha and beta functions
    phi = Q_10^((T_actual - T_base)/10); % (dimless) temperature scaling factor
    alpha_n = @(Vm) phi * 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10)); % (1/ms)
    beta_n = @(Vm) phi * 0.125*exp(-(Vm + 65)/80); % (1/ms)
    alpha_m = @(Vm) phi * 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10)); % (1/ms)
    beta_m = @(Vm) phi * 4*exp(-(Vm + 65)/18); % (1/ms)
    alpha_h = @(Vm) phi * 0.07*exp(-(Vm + 65)/20); % (1/ms)
    beta_h = @(Vm) phi * 1/(1 + exp(-(Vm + 35)/10)); % (1/ms)

    % Defining the initial data
    Vm = V_m0 * ones(1, m);
    N = N_0 * ones(1, m);
    M = M_0 * ones(1, m);
    H = H_0 * ones(1, m);
    
    % Defining the matrices to collect data at every time step
    Vm_all(1,:) = Vm;
    N_all(1,:) = N;
    M_all(1,:) = M; 
    H_all(1,:) = H;
    
    % Defining the initial matrix A and vectors g_1 and g_2
    A = zeros(m);
    g_1 = zeros(m, 1);
    g_2 = zeros(m, 1);

    % Starting the time loop
    for j = 1:(n-1)
        
        % setting newN, newM, newH vectors
        for i = 1:m
            newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
            newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
            newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
        end
        
        % Edit the next Vm, N, M and H (redefining Vm, N, M and H vectors)
        N = newN;
        M = newM;
        H = newH;
        
        % Setting the boundary conditions
        A(1, 1) = 1;
        A(1, 2) = -1;
        A(m, m-1) = -1;
        A(m, m) = 1;
        
        % Starting the space loop
        for i = 2:m-1
    
            gamma1 = -dt*a/(2*R_i*C_m*dx^2);
            gamma2 = 1 + dt*a/(R_i*C_m*dx^2) + dt*G_K*(N(i)^4)/C_m + dt*G_Na*(M(i)^3)*H(i)/C_m + dt*G_L/C_m;
            gamma3 = -dt*a/(2*R_i*C_m*dx^2); 
            gamma4 = 1; 
            gamma5 = dt*G_K*N(i)^4*E_K/C_m + dt*G_Na*(M(i)^3)*H(i)*E_Na/C_m + dt*G_L*E_L/C_m;
    
            % adding stimulus temporally and spatially:
            if (j*dt >= S_T0 && j*dt <= S_T1) && (i*dx >= S_P0 && i*dx <= S_P1)
                gamma2 = 1 + dt*a/(R_i*C_m*dx^2) + dt*G_K*(N(i)^4)/C_m + dt*(G_Na*(M(i)^3)*H(i) + S_v)/C_m + dt*G_L/C_m;
                gamma5 = dt*G_K*(N(i)^4)*E_K/C_m + dt*(G_Na*(M(i)^3)*H(i) + S_v)*E_Na/C_m + dt*G_L*E_L/C_m;
            end
            
            A(i, i-1) = gamma3;
            A(i, i) = gamma2;
            A(i, i+1) = gamma1;
            g_1(i) = gamma4*Vm(i); 
            g_2(i) = gamma5;
            
        end
        
        % Setting newVm
        newVm = transpose(A\(g_1+g_2));
        Vm = newVm;
        
        % Adding the newly defined vectors to the 'all' matrices
        Vm_all(j+1,:) = Vm;
        N_all(j+1,:) = N;
        M_all(j+1,:) = M;
        H_all(j+1,:) = H;
    
        j % displaying the current time iteration (nice way to see how long simulation is)
    
    end

    % Saving all the data (can be saved outside of the function in main)
    hh_solution_data.Vm_all = Vm_all;
    hh_solution_data.N_all = N_all;
    hh_solution_data.M_all = M_all;
    hh_solution_data.H_all = H_all;
    hh_solution_data.mesh_params = mesh_params;
    hh_solution_data.material_params = material_params;
    
    % Saving all the data defined in this function (automatically saved)
    save('hh_simulation.mat');

end
% SC model v3 scheme
% Kevin Roberts
% April 2025

% This function will solve for Vm, n, m and h for the SC model given mesh
% and material parameters

function sc_solution_data = sc_function_v3(mesh_params, material_params)
    
    % Grabing and defining all the inputed mesh and material parameters
    a = mesh_params.a; 
    a_my = mesh_params.a_my;
    dx = mesh_params.dx;
    dt = mesh_params.dt; 
    L = mesh_params.L;
    T = mesh_params.T;
    N_n = mesh_params.N_n;
    N_s = mesh_params.N_s;
    m = mesh_params.m;
    n = mesh_params.n;
    R_i = material_params.R_i;
    R_m = material_params.R_m;
    C_m = material_params.C_m;  
    R_my = material_params.R_my; 
    C_my = material_params.C_my; 
    G_K = material_params.G_K; 
    G_Na = material_params.G_Na; 
    G_L = material_params.G_L; 
    E_K = material_params.E_K;
    E_Na = material_params.E_Na;
    E_L = material_params.E_L; 
    E_rest = material_params.E_rest;
    S_v = material_params.S_v;
    S_T0 = material_params.S_T0;
    S_T1 = material_params.S_T1;
    S_P0 = material_params.S_P0;
    S_P1 = material_params.S_P1;
    T_base = material_params.T_base;
    T_actual = material_params.T_actual;
    Q_10_Na = material_params.Q_10_Na;
    Q_10_K = material_params.Q_10_K; 
    V_m0 = material_params.V_m0; 
    V_my0 = material_params.V_my0; 
    N_0 = material_params.N_0; 
    M_0 = material_params.M_0; 
    H_0 = material_params.H_0; 

    % Defining rho and w_1 constants
    rho = dt/dx^2;
    w_1 = a^2/(C_my*a_my*R_i);
    
    % Defining the stimulus function. Note that in the S function: ii is the space index and tt is the time index
    S = @(ii, tt) S_v * ((abs(tt * dt - S_T0) <= 1e-10 | tt * dt > S_T0) & ...
                        (tt * dt < S_T1 | abs(tt * dt - S_T1) <= 1e-10) & ...
                        (abs(ii * dx - S_P0) <= 1e-10 | ii * dx > S_P0) & ...
                        (ii * dx < S_P1 | abs(ii * dx - S_P1) <= 1e-10));
    
    % Defining alpha/beta functions
    phi_Na = Q_10_Na^((T_actual - T_base)/10); % (#) temperature scaling factor for Na current
    phi_K = Q_10_K^((T_actual - T_base)/10); % (#) temperature scaling factor for K current
    alpha_n = @(Vm) phi_K * 0.01*(Vm + 55)/(1 - exp(-(Vm + 55)/10));
    beta_n = @(Vm) phi_K * 0.125*exp(-(Vm + 65)/80);
    alpha_m = @(Vm) phi_Na * 0.1*(Vm + 40)/(1 - exp(-(Vm + 40)/10));
    beta_m = @(Vm) phi_Na * 4*exp(-(Vm + 65)/18);
    alpha_h = @(Vm) phi_Na * 0.07*exp(-(Vm + 65)/20);
    beta_h = @(Vm) phi_Na * 1/(1 + exp(-(Vm + 35)/10));
    
    % Defining the b_1(x_i) function
    B_1 = (a/(2*R_i*C_m))*(1 + C_m*a/(C_my*a_my)); % Internodal region
    B_2 = a/(2*R_i*C_m); % Nodal region
    B_3 = (B_1 + B_2)/2; % End point (may not be used)
    b_1 = @(ii) (mod(ii - 1, N_s) > N_n).*B_1 + ... % Internodal region
               (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*B_2 + ... % Nodal region
               ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*B_3; % End point      
    
    % Defining the c_1(x_i) function
    C_1 = -1/(R_m*C_m); % Internodal region
    C_2 = @(n, m, h, ii, tt) -1/C_m*(G_K*n^4 + (G_Na*m^3*h + S(ii, tt)) + G_L); % Nodal region
    C_3 = @(n, m, h, ii, tt) (C_1 + C_2(n, m, h, ii, tt))/2; % End point
    c_1 = @(n, m, h, ii, tt) (mod(ii - 1, N_s) > N_n).*C_1 + ... % Internodal region
               (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*C_2(n, m, h, ii, tt) + ... % Nodal region
               ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*C_3(n, m, h, ii, tt); % End point        
           
    % Defining the f_1(x_i) function
    F_1 = @(Vmy) (1/(R_m*C_m) - 1/(C_my*R_my))*Vmy + E_rest/(R_m*C_m); % Internodal region
    F_2 = @(n, m, h, ii, tt) 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(ii, tt))*E_Na + G_L*E_L); % Nodal region
    F_3 = @(Vmy, n, m, h, ii, tt) (F_1(Vmy) + F_2(n, m, h, ii, tt))/2; % End point
    f_1 = @(Vmy, n, m, h, ii, tt) (mod(ii - 1, N_s) > N_n).*F_1(Vmy) + ... % Internodal region
               (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*F_2(n, m, h, ii, tt) + ... % Nodal region
               ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*F_3(Vmy, n, m, h, ii, tt); % End point        
    
    % Defining the initial vectors
    Vm = V_m0 * ones(1, m);
    Vmy = zeros(1, m);
    N = zeros(1, m);
    M = zeros(1, m);
    H = zeros(1, m);
    
    % Populating n, m h, and Vmy with appropriate initial conditions along
    % the myelinated axon
    N(1) = N_0;
    M(1) = M_0;
    H(1) = H_0;
    for i = 2:m-1
        
        seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
        myelin_start = (seg - 1)*(N_s) + N_n; % start of internodal region in this segment
        myelin_end = seg*(N_s); % end of internodal region in this segment
        seg_start = (seg - 1)*(N_s); % index of the start of the segment
    
        % Internodal region
        if (i > myelin_start + 1) && (i < myelin_end + 1)
            Vmy(i) = V_my0;
            N(i) = 0;
            M(i) = 0;
            H(i) = 0;
        % Nodal region
        elseif (i > seg_start + 1) && (i < myelin_start + 1)
            Vmy(i) = 0;
            N(i) = N_0;
            M(i) = M_0;
            H(i) = H_0;
        % End points
        else
            Vmy(i) = 0;
            N(i) = N_0;
            M(i) = M_0;
            H(i) = H_0;
        end 
    end
    Vmy(m) = 0;
    
    % Defining the matrices to collect data at every time step
    Vm_all(1,:) = Vm;
    Vmy_all(1,:) = Vmy;
    Vm_minus_Vmy(1,:) = Vm - Vmy;
    N_all(1,:) = N;
    M_all(1,:) = M; 
    H_all(1,:) = H;
    
    % Running the time loop
    %%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:(n-1)
        
        % Updating the probability gate functions n, m and h
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:m-1
            seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
            myelin_start = (seg - 1)*(N_s) + N_n; % start of internodal region in this segment
            myelin_end = seg*(N_s); % end of internodal region in this segment
            seg_start = (seg - 1)*(N_s); % index of the start of the segment
    
            if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
                newN(i) = 0;
                newM(i) = 0;
                newH(i) = 0;
            elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
                newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
                newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
                newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
            else % End point
                newN(i) = 1/(1 + dt*alpha_n(Vm(i)) + dt*beta_n(Vm(i))) * (N(i) + dt*alpha_n(Vm(i)));
                newM(i) = 1/(1 + dt*alpha_m(Vm(i)) + dt*beta_m(Vm(i))) * (M(i) + dt*alpha_m(Vm(i)));
                newH(i) = 1/(1 + dt*alpha_h(Vm(i)) + dt*beta_h(Vm(i))) * (H(i) + dt*alpha_h(Vm(i)));
            end
        end
        newN(m) = 0;
        newM(m) = 0;
        newH(m) = 0;
        

        % Updating Vmy
        %%%%%%%%%%%%%%
        newVmy(1) = 0;
        for i = 2:m-1
            seg = floor((i - 1)/(N_s)) + 1; % axon segment number based on index i
            myelin_start = (seg - 1)*(N_s) + N_n; % start of internodal region in this segment
            myelin_end = seg*(N_s); % end of internodal region in this segment
    
            if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
                
                eta1 = 1/(1 + dt/(C_my*R_my));
                eta2 = rho * w_1/2 * 1/(1 + dt/(C_my*R_my));
                eta3 = -rho * w_1 * 1/(1 + dt/(C_my*R_my));
                eta4 = rho * w_1/2 * 1/(1 + dt/(C_my*R_my));
            
            else % End point and Nodal Regions
                eta1 = 0;
                eta2 = 0;
                eta3 = 0;
                eta4 = 0;
            end
            newVmy(i) = eta1*Vmy(i) + eta2*Vm(i-1) + eta3*Vm(i) + eta4*Vm(i+1);
        
        end
        newVmy(m) = 0;
        

        % Updating Vm

        % Defining the A matrix and g_1 and g_2 vectors
        A = zeros(m, m);
        g_1 = zeros(m, 1);
        g_2 = zeros(m, 1);
        
        % Setting the boundary conditions
        A(1, 1) = 1;
        A(1, 2) = -1;
        A(m, m-1) = -1;
        A(m, m) = 1;
    
        % Starting the space loop
        for i = 2:(m-1)
    
            gamma1 = -rho*b_1(i - 1/2);
            gamma2 = 1 - dt*c_1(newN(i), newM(i), newH(i), i, j) + rho*(b_1(i + 1/2) + b_1(i - 1/2));
            gamma3 = -rho*b_1(i + 1/2);
            gamma4 = 1;
            gamma5 = dt * f_1(newVmy(i), newN(i), newM(i), newH(i), i, j);
    
            A(i, i-1) = gamma1;
            A(i, i) = gamma2;
            A(i, i+1) = gamma3;
            g_1(i) = gamma4*Vm(i);
            g_2(i) = gamma5;
               
        end

        j % displaying the current time iteration (nice way to see how long simulation is)

        % Solving for V_m^{j+1}
        %%%%%%%%%%%%%%%%%%%%%%%
        newVm = transpose(A\(g_1+g_2));
    
        % Updating Vm, Vmy, N, M, and H
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Vm = newVm;
        Vmy = newVmy;
        N = newN;
        M = newM;
        H = newH;
        
        % Adding the updated vectors to the 'all' data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Vm_all(j+1,:) = Vm;
        Vmy_all(j+1,:) = Vmy;
        Vm_minus_Vmy(j+1,:) = Vm - Vmy;
        N_all(j+1,:) = N;
        M_all(j+1,:) = M;
        H_all(j+1,:) = H;
    
    end

    % Saving all the data (can be saved outside of the function in main)
    sc_solution_data.Vm_all = Vm_all;
    sc_solution_data.Vmy_all = Vmy_all;
    sc_solution_data.N_all = N_all;
    sc_solution_data.M_all = M_all;
    sc_solution_data.H_all = H_all;
    sc_solution_data.mesh_params = mesh_params;
    sc_solution_data.material_params = material_params;
    
    % Saving all the data defined in this function (automatically saved)
    save('sc_simulation_v3.mat');

end
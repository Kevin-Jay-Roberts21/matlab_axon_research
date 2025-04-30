% DC model v3 scheme
% Kevin Roberts
% April 2025

% This function will solve for Vm, n, m and h for the DC model given mesh
% and material parameters

function dc_solution_data = dc_function_v3(mesh_params, material_params)
    
    % Grabing and defining all the inputed mesh and material parameters
    dx = mesh_params.dx;
    dt = mesh_params.dt; 
    L_my = mesh_params.L_my;
    L_n = mesh_params.L_n;
    L_pn = mesh_params.L_pn;
    d_pa = mesh_params.d_pa;
    d_pn = mesh_params.d_pn;
    L_s = mesh_params.L_s;
    n_s = mesh_params.n_s;
    L = mesh_params.L;
    T = mesh_params.T;
    N_n = mesh_params.N_n;
    N_my = mesh_params.N_my;
    N_s = mesh_params.N_s;
    m = mesh_params.m;
    n = mesh_params.n;
    a = material_params.a; 
    a_my = material_params.a_my;
    R_i = material_params.R_i;
    R_m = material_params.R_m;
    C_m = material_params.C_m; 
    r_pa = material_params.r_pa; 
    r_pn = material_params.r_pn; 
    R_my = material_params.R_my; 
    C_my = material_params.C_my; 
    R_pa = material_params.R_pa; 
    R_pn = material_params.R_pn; 
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

    % Defining rho, w1, w2 and w3 constants
    rho = dt/dx^2; % creating the courant number
    w1 = a^2/(C_my*a_my*R_i);
    w2 = d_pa*(2*a + d_pa)/(C_my*a_my*R_pa);
    w3 = r_pa/(r_pn*L_pn);

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
           
    % defining the f_2(x_i) function
    F_1 = @(Vmy_i) (1/(R_m*C_m) - 1/(C_my*R_my))*Vmy_i + E_rest/(R_m*C_m);
    F_4 = @(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1) w2/(2*dx^2)*Vmy_i_minus_1 - w2/dx^2*Vmy_i + w2/(2*dx^2)*Vmy_i_plus_1; % Internodal region
    F_2 = @(n, m, h, ii, tt) 1/C_m*(G_K*n^4*E_K + (G_Na*m^3*h + S(ii, tt))*E_Na + G_L*E_L); % Nodal region
    F_5 = @(Vmy_i, n, m, h, ii, tt) (F_1(Vmy_i) + F_2(n, m, h, ii, tt))/2; % End point
     
    f_2 = @(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1, n, m, h, ii, tt) (mod(ii - 1, N_s) > N_n).*(F_1(Vmy_i) + F_4(Vmy_i_minus_1, Vmy_i, Vmy_i_plus_1)) + ... % Internodal region
               (mod(ii - 1, N_s) < N_n & mod(ii - 1, N_s) ~= 0).*F_2(n, m, h, ii, tt) + ... % Nodal region
               ((mod(ii - 1, N_s) == N_n) | (mod(ii - 1, N_s) == 0)).*F_5(Vmy_i, n, m, h, ii, tt); % End point  
    
    % Defining the initial vectors
    Vm = V_m0 * ones(1, m);
    Vmy = V_my0 * ones(1, m);
    N = zeros(1, m);
    M = zeros(1, m);
    H = zeros(1, m);
    
    % Populating n, m h, and Vmy with appropriate initial conditions along
    % the myelinated axon
    Vmy(1) = 0;
    N(1) = N_0;
    M(1) = M_0;
    H(1) = H_0;
    for i = 2:m-1 % because we want to keep the end points (1 and M) for Vmy, n, m, h at 0
        
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
        elseif (i == seg_start + 1) % Right end point (x_R)
            Vmy(i) = V_my0/(1 + w3*dx);
            N(i) = N_0;
            M(i) = M_0;
            H(i) = H_0;
        elseif (i == myelin_start + 1) % Left end point (x_L)
            Vmy(i) = V_my0/(1 + w3*dx);
            N(i) = N_0;
            M(i) = M_0;
            H(i) = H_0;
        end 
    end
    Vmy(m) = V_my0/(1 + w3*dx);
    
    % Defining the matrices to collect data at every time step
    Vm_all(1,:) = Vm;
    Vmy_all(1,:) = Vmy;
    Vm_minus_Vmy(1,:) = Vm - Vmy;
    N_all(1,:) = N;
    M_all(1,:) = M; 
    H_all(1,:) = H;

    % Running the time loop
    for j = 1:(n-1)
        
        % updating the probability gate functions n, m and h
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
        
        % Defining the A_2 matrix as well as e_2 and f_2
        A_2 = zeros(m, m);
        e_1 = zeros(m, 1);
        e_2 = zeros(m, 1);
    
        % Using the boundary conditions to define the top and bottom row of A_2
        A_2(1, 1) = 1;
        A_2(1, 2) = 0;
        A_2(m, m-1) = -1;
        A_2(m, m) = (1 + w3*dx);
    
        % Starting the space loop
        for i = 2:(m-1)
            % updating the A_2 matrix each row is determined by the
            % internodal, nodal, and end points (left or right end points)
            seg = floor((i - 1)/N_s) + 1; % axon segment number based on index i
            myelin_start = (seg - 1)*N_s + N_n; % start of internodal region in this segment
            myelin_end = seg*N_s; % end of internodal region in this segment
            seg_start = (seg - 1)*N_s; % index of the start of the segment
    
            if (i > myelin_start + 1) && (i < myelin_end + 1) % Internodal region
                eta1 = -rho*w2/2;
                eta2 = 1 + dt/(R_my*C_my) + rho*w2;
                eta3 = -rho*w2/2; 
                eta4 = 1; 
                eta5 = rho*w1/2*Vm(i-1) - rho*w1*Vm(i) + rho*w1/2*Vm(i+1);
            elseif (i > seg_start + 1) && (i < myelin_start + 1) % Nodal region
                eta1 = 0;
                eta2 = 1;
                eta3 = 0; 
                eta4 = 0; 
                eta5 = 0;
            elseif (i == seg_start + 1) % Right end point (x_R)
                eta1 = -1;
                eta2 = (1 + dx*w3);
                eta3 = 0; 
                eta4 = 0; 
                eta5 = 0;
            elseif (i == myelin_start + 1) % Left end point (x_L)
                eta1 = 0;
                eta2 = (1 + dx*w3);
                eta3 = -1;
                eta4 = 0; 
                eta5 = 0;
            end
    
            A_2(i, i-1) = eta1;
            A_2(i, i) = eta2;
            A_2(i, i+1) = eta3;
            e_1(i) = eta4*Vmy(i);
            e_2(i) = eta5;
    
        end
        newVmy = transpose(A_2\(e_1+e_2));
        

        % Updating Vm 
    
        % Defining the A_1 matrix
        A_1 = zeros(m, m);
        g_1 = zeros(m, 1);
        g_2 = zeros(m, 1);
    
        % using the boundary conditions to define the top and bottom row of A
        A_1(1, 1) = 1;
        A_1(1, 2) = -1;
        A_1(m, m-1) = -1;
        A_1(m, m) = 1;
        
        % starting the space loop
        for i = 2:(m-1)
    
            gamma1 = -rho*b_1(i - 1/2);
            gamma2 = 1 - dt*c_1(newN(i), newM(i), newH(i), i, j) + rho*(b_1(i + 1/2) + b_1(i - 1/2));
            gamma3 = -rho*b_1(i + 1/2);
            gamma4 = 1;
            gamma5 = dt * f_2(newVmy(i-1), newVmy(i), newVmy(i+1), newN(i), newM(i), newH(i), i, j);
    
            A_1(i, i-1) = gamma1;
            A_1(i, i) = gamma2;
            A_1(i, i+1) = gamma3;
            g_1(i) = gamma4*Vm(i);
            g_2(i) = gamma5;
            
        end
        
        j % displaying the current time iteration (nice way to see how long simulation is)
    
        % Solving for V_m^{j+1}
        newVm = transpose(A_1\(g_1+g_2));
    
        % Updating Vmy and Vm and adding the data to the _all matrices
        Vm = newVm;
        Vmy = newVmy;
        N = newN;
        M = newM;
        H = newH;
        
        % Adding the updated vectors to the '_all' data
        Vm_all(j+1,:) = Vm;
        Vmy_all(j+1,:) = Vmy;
        Vm_minus_Vmy(j+1,:) = Vm - Vmy;
        N_all(j+1,:) = N;
        M_all(j+1,:) = M;
        H_all(j+1,:) = H;
    
    end


    % Saving all the data (can be saved outside of the function in main)
    dc_solution_data.Vm_all = Vm_all;
    dc_solution_data.Vmy_all = Vmy_all;
    dc_solution_data.N_all = N_all;
    dc_solution_data.M_all = M_all;
    dc_solution_data.H_all = H_all;
    dc_solution_data.mesh_params = mesh_params;
    dc_solution_data.material_params = material_params;
    
    % Saving all the data defined in this function (automatically saved)
    save('dc_simulation_v3.mat');

end
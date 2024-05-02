function [Linf_U, Linf_N, Linf_M, Linf_H, L1_U, L1_N, L1_M, L1_H, L2_U, L2_N, L2_M, L2_H] = norm_calculator(U1, N1, M1, H1, U2, N2, M2, H2, num_of_rows)

    
    % compute the Linf norm difference of -1s and -2s
    % compute the L1 norm difference of -1s and -2s
    % compute the L2 norm difference of -1s and -2s
    U_max_numerator = 0;
    U_max_denominator = 0;
    N_max_numerator = 0;
    N_max_denominator = 0;
    M_max_numerator = 0; 
    M_max_denominator = 0; 
    H_max_numerator = 0;
    H_max_denominator = 0;

    U_L1_norm_numerator = 0;
    U_L1_norm_denominator = 0;
    N_L1_norm_numerator = 0;
    N_L1_norm_denominator = 0;
    M_L1_norm_numerator = 0;
    M_L1_norm_denominator = 0;
    H_L1_norm_numerator = 0;
    H_L1_norm_denominator = 0;

    U_L2_norm_numerator_sum = 0;
    U_L2_norm_denominator_sum = 0;
    N_L2_norm_numerator_sum = 0;
    N_L2_norm_denominator_sum = 0;
    M_L2_norm_numerator_sum = 0;
    M_L2_norm_denominator_sum = 0;
    H_L2_norm_numerator_sum = 0;
    H_L2_norm_denominator_sum = 0;


    for i = 1:num_of_rows

        % finding the L1 norm
        U_L1_norm_numerator = U_L1_norm_numerator + abs(U1(i, 1) - U2(i, 1));
        U_L1_norm_denominator = U_L1_norm_denominator + abs(U2(i, 1));
        N_L1_norm_numerator = N_L1_norm_numerator + abs(N1(i, 1) - N2(i, 1));
        N_L1_norm_denominator = N_L1_norm_denominator + abs(N2(i, 1));
        M_L1_norm_numerator = M_L1_norm_numerator + abs(M1(i, 1) - M2(i, 1));
        M_L1_norm_denominator = M_L1_norm_denominator + abs(M2(i, 1));
        H_L1_norm_numerator = H_L1_norm_numerator + abs(H1(i, 1) - H2(i, 1));
        H_L1_norm_denominator = H_L1_norm_denominator + abs(H2(i, 1));

        % finding the L2 norm
        U_L2_norm_numerator_sum = U_L2_norm_numerator_sum + (abs(U1(i, 1) - U2(i, 1)))^2;
        U_L2_norm_denominator_sum = U_L2_norm_denominator_sum + (abs(U2(i, 1)))^2;
        N_L2_norm_numerator_sum = N_L2_norm_numerator_sum + (abs(N1(i, 1) - N2(i, 1)))^2;
        N_L2_norm_denominator_sum = N_L2_norm_denominator_sum + (abs(N2(i, 1)))^2;
        M_L2_norm_numerator_sum = M_L2_norm_numerator_sum + (abs(M1(i, 1) - M2(i, 1)))^2;
        M_L2_norm_denominator_sum = M_L2_norm_denominator_sum + (abs(M2(i, 1)))^2;
        H_L2_norm_numerator_sum = H_L2_norm_numerator_sum + (abs(H1(i, 1) - H2(i, 1)))^2;
        H_L2_norm_denominator_sum = H_L2_norm_denominator_sum + (abs(H2(i, 1)))^2;    

        % finding the Linf norm
        new_umax_numerator = abs(U1(i, 1) - U2(i, 1));
        new_nmax_numerator = abs(N1(i, 1) - N2(i, 1));
        new_mmax_numerator = abs(M1(i, 1) - M2(i, 1));
        new_hmax_numerator = abs(H1(i, 1) - H2(i, 1));
        new_umax_denominator = abs(U2(i, 1));
        new_nmax_denominator = abs(N2(i, 1));
        new_mmax_denominator = abs(M2(i, 1));
        new_hmax_denominator = abs(H2(i, 1));

        % find Linf norm
        if new_umax_numerator > U_max_numerator
            U_max_numerator = new_umax_numerator;
        end
        if new_nmax_numerator > N_max_numerator
            N_max_numerator = new_nmax_numerator;
        end
        if new_mmax_numerator > M_max_numerator
            M_max_numerator = new_mmax_numerator;
        end
        if new_hmax_numerator > H_max_numerator
            H_max_numerator = new_hmax_denominator;
        end

        if new_umax_denominator > U_max_denominator
            U_max_denominator = new_umax_denominator;
        end
        if new_nmax_denominator > N_max_denominator
            N_max_denominator = new_nmax_denominator;
        end
        if new_mmax_denominator > M_max_denominator
            M_max_denominator = new_mmax_denominator;
        end
        if new_hmax_denominator > H_max_denominator
            H_max_denominator = new_hmax_denominator;
        end

    end

    % Linf norms
    Linf_U = round(vpa(U_max_numerator/U_max_denominator), 5);
    Linf_N = round(vpa(N_max_numerator/N_max_denominator), 5);
    Linf_M = round(vpa(M_max_numerator/M_max_denominator), 5);
    Linf_H = round(vpa(H_max_numerator/H_max_denominator), 5);

    % L1 norms
    L1_U = round(vpa(U_L1_norm_numerator/U_L1_norm_denominator), 5);
    L1_N = round(vpa(N_L1_norm_numerator/N_L1_norm_denominator), 5);
    L1_M = round(vpa(M_L1_norm_numerator/M_L1_norm_denominator), 5);
    L1_H = round(vpa(H_L1_norm_numerator/H_L1_norm_denominator), 5);

    % L2 norms
    L2_U = round(vpa(sqrt(U_L2_norm_numerator_sum)/sqrt(U_L2_norm_denominator_sum)), 5);
    L2_N = round(vpa(sqrt(N_L2_norm_numerator_sum)/sqrt(N_L2_norm_denominator_sum)), 5);
    L2_M = round(vpa(sqrt(M_L2_norm_numerator_sum)/sqrt(M_L2_norm_denominator_sum)), 5);
    L2_H = round(vpa(sqrt(H_L2_norm_numerator_sum)/sqrt(H_L2_norm_denominator_sum)), 5);
    
end


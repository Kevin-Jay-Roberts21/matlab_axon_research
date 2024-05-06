function [Linf_U, Linf_N, Linf_M, Linf_H, L1_U, L1_N, L1_M, L1_H, L2_U, L2_N, L2_M, L2_H] = norm_calculator(U1, N1, M1, H1, U2, N2, M2, H2, num_of_rows)

    % compute the Linf norm difference of -1s and -2s
    % compute the L1 norm difference of -1s and -2s
    % compute the L2 norm difference of -1s and -2s

    % Linf norms
    Linf_U = round(max(abs(U1-U2))/max(abs(U2)), 5);
    Linf_N = round(max(abs(N1-N2))/max(abs(N2)), 5);
    Linf_M = round(max(abs(M1-M2))/max(abs(M2)), 5);
    Linf_H = round(max(abs(H1-H2))/max(abs(H2)), 5);

    % L1 norms
    L1_U = round(sum(abs(U1-U2))/sum(abs(U2)), 5);
    L1_N = round(sum(abs(N1-N2))/sum(abs(N2)), 5);
    L1_M = round(sum(abs(M1-M2))/sum(abs(M2)), 5);
    L1_H = round(sum(abs(H1-H2))/sum(abs(H2)), 5);
    
    % L2 norms
    L2_U = round(sqrt(sum(abs(U1-U2).^2)/sum(abs(U2).^2)), 5);
    L2_N = round(sqrt(sum(abs(N1-N2).^2)/sum(abs(N2).^2)), 5);
    L2_M = round(sqrt(sum(abs(M1-M2).^2)/sum(abs(M2).^2)), 5);
    L2_H = round(sqrt(sum(abs(H1-H2).^2)/sum(abs(H2).^2)), 5);
    
end


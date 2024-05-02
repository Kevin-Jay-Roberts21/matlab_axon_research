clear all
close all
clc

bigUall1 = load('projects/axon_simulations/data/Uall1.mat').Uall;
bigNall1 = load('projects/axon_simulations/data/Nall1.mat').Nall;
bigMall1 = load('projects/axon_simulations/data/Mall1.mat').Mall;
bigHall1 = load('projects/axon_simulations/data/Hall1.mat').Hall;
Uall1_num_of_rows = size(bigUall1, 1);
Uall1_num_of_cols = size(bigUall1, 2);

bigUall2 = load('projects/axon_simulations/data/Uall2.mat').Uall;
bigNall2 = load('projects/axon_simulations/data/Nall2.mat').Nall;
bigMall2 = load('projects/axon_simulations/data/Mall2.mat').Mall;
bigHall2 = load('projects/axon_simulations/data/Hall2.mat').Hall;
Uall2_num_of_rows = size(bigUall2, 1);
Uall2_num_of_cols = size(bigUall2, 2);

bigUall3 = load('projects/axon_simulations/data/Uall3.mat').Uall;
bigNall3 = load('projects/axon_simulations/data/Nall3.mat').Nall;
bigMall3 = load('projects/axon_simulations/data/Mall3.mat').Mall;
bigHall3 = load('projects/axon_simulations/data/Hall3.mat').Hall;
Uall3_num_of_rows = size(bigUall3, 1);
Uall3_num_of_cols = size(bigUall3, 2);

U1 = zeros(Uall1_num_of_cols, 1);
N1 = zeros(Uall1_num_of_cols, 1);
M1 = zeros(Uall1_num_of_cols, 1);
H1 = zeros(Uall1_num_of_cols, 1);
U2 = zeros(Uall1_num_of_cols, 1);
N2 = zeros(Uall1_num_of_cols, 1);
M2 = zeros(Uall1_num_of_cols, 1);
H2 = zeros(Uall1_num_of_cols, 1);
U3 = zeros(Uall1_num_of_cols, 1);
N3 = zeros(Uall1_num_of_cols, 1);
M3 = zeros(Uall1_num_of_cols, 1);
H3 = zeros(Uall1_num_of_cols, 1);


% redefining Uall1, Uall2 and Uall3 to be vectors of the same length
for i = 1:Uall1_num_of_rows
    U1(i, 1) = bigUall1(i, 1);
    N1(i, 1) = bigNall1(i, 1);
    M1(i, 1) = bigMall1(i, 1);
    H1(i, 1) = bigHall1(i, 1);
    
    U2(i, 1) = bigUall2((2*i)-1, 1);
    N2(i, 1) = bigNall2((2*i)-1, 1);
    M2(i, 1) = bigMall2((2*i)-1, 1);
    H2(i, 1) = bigHall2((2*i)-1, 1);
    
    U3(i, 1) = bigUall3(10*(i-1)+1, 1);
    N3(i, 1) = bigNall3(10*(i-1)+1, 1);
    M3(i, 1) = bigMall3(10*(i-1)+1, 1);
    H3(i, 1) = bigHall3(10*(i-1)+1, 1);
end

% computing h = 0.01 vs h = 0.005 norms
[Linf_U, Linf_N, Linf_M, Linf_H, L1_U, L1_N, L1_M, L1_H, L2_U, L2_N, L2_M, L2_H] = norm_calculator(U1, N1, M1, H1, U2, N2, M2, H2, Uall1_num_of_rows)

% computing H = 0.005 vs h = 0.001 norms
[Linf_U, Linf_N, Linf_M, Linf_H, L1_U, L1_N, L1_M, L1_H, L2_U, L2_N, L2_M, L2_H] = norm_calculator(U2, N2, M2, H2, U3, N3, M3, H3, Uall1_num_of_rows)


% plotting h = 0.01, h= 0.005, and h = 0.001 next to each other
figure(1)
plot_different_time_steps(U1, U2, U3, 0, 35, 0.01)

% plotting from time t = 6 milliseconds to t = 7 milliseconds
figure(2)
plot_different_time_steps(U1, U2, U3, 6, 7, 0.01)

% plotting from time t = 5 milliseconds to t = 10 milliseconds
figure(3)
plot_different_time_steps(U1, U2, U3, 5, 10, 0.01)

clear all
close all
clc

Uall1 = load('Uall1.mat').Uall;
Nall1 = load('Nall1.mat').Nall;
Mall1 = load('Mall1.mat').Mall;
Hall1 = load('Hall1.mat').Hall;
Uall1_num_of_rows = size(Uall1, 1);
Uall1_num_of_cols = size(Uall1, 2);

Uall2 = load('Uall2.mat').Uall;
Nall2 = load('Nall2.mat').Nall;
Mall2 = load('Mall2.mat').Mall;
Hall2 = load('Hall2.mat').Hall;
Uall2_num_of_rows = size(Uall2, 1);
Uall2_num_of_cols = size(Uall2, 2);

Uall3 = load('Uall3.mat').Uall;
Nall3 = load('Nall3.mat').Nall;
Mall3 = load('Mall3.mat').Mall;
Hall3 = load('Hall3.mat').Hall;
Uall3_num_of_rows = size(Uall3, 1);
Uall3_num_of_cols = size(Uall3, 2);


% find the difference at time t_0 = 0.01, t_1 = 0.02, ..., t_3500 = 35
% Modify Uall2, Nall2, Mall2 Hall2, and all the others to have the same
% number of rows as Uall1, Nall1, Mall1 and Hall1

% compute the Linf norm difference of -all1s and -all2s
% compute the L1 norm difference of -all1s and -all2s
% compute the L2 norm difference of -all1s and -all2s
Uall_max_numerator = 0;
Nall_max_numerator = 0;
Mall_max_numerator = 0; 
Hall_max_numerator = 0;
Uall_max_denominator = 0;
Nall_max_denominator = 0;
Mall_max_denominator = 0; 
Hall_max_denominator = 0;
Uall_L1_norm_numerator = 0;
Uall_L1_norm_denominator = 0;
Uall_L2_norm_numerator_sum = 0;
Uall_L2_norm_denominator_sum = 0;
Nall_L1_norm_numerator = 0;
Nall_L1_norm_denominator = 0;
Nall_L2_norm_numerator_sum = 0;
Nall_L2_norm_denominator_sum = 0;
Mall_L1_norm_numerator = 0;
Mall_L1_norm_denominator = 0;
Mall_L2_norm_numerator_sum = 0;
Mall_L2_norm_denominator_sum = 0;
Hall_L1_norm_numerator = 0;
Hall_L1_norm_denominator = 0;
Hall_L2_norm_numerator_sum = 0;
Hall_L2_norm_denominator_sum = 0;

for i = 1:Uall1_num_of_rows
    
    % finding the L1 norm
    Uall_L1_norm_numerator = Uall_L1_norm_numerator + abs(Uall1(i, 1) - Uall2((2*i)-1, 1));
    Uall_L1_norm_denominator = Uall_L1_norm_denominator + abs(Uall2((2*i)-1, 1));
    Nall_L1_norm_numerator = Nall_L1_norm_numerator + abs(Nall1(i, 1) - Nall2((2*i)-1, 1));
    Nall_L1_norm_denominator = Nall_L1_norm_denominator + abs(Nall2((2*i)-1, 1));
    Mall_L1_norm_numerator = Mall_L1_norm_numerator + abs(Mall1(i, 1) - Mall2((2*i)-1, 1));
    Mall_L1_norm_denominator = Mall_L1_norm_denominator + abs(Mall2((2*i)-1, 1));
    Hall_L1_norm_numerator = Hall_L1_norm_numerator + abs(Hall1(i, 1) - Hall2((2*i)-1, 1));
    Hall_L1_norm_denominator = Hall_L1_norm_denominator + abs(Hall2((2*i)-1, 1));
    
    % finding the L2 norm
    Uall_L2_norm_numerator_sum = Uall_L2_norm_numerator_sum + (abs(Uall1(i, 1) - Uall2((2*i)-1, 1)))^2;
    Uall_L2_norm_denominator_sum = Uall_L2_norm_denominator_sum + (abs(Uall2((2*i)-1, 1)))^2;
    Nall_L2_norm_numerator_sum = Nall_L2_norm_numerator_sum + (abs(Nall1(i, 1) - Nall2((2*i)-1, 1)))^2;
    Nall_L2_norm_denominator_sum = Nall_L2_norm_denominator_sum + (abs(Nall2((2*i)-1, 1)))^2;
    Mall_L2_norm_numerator_sum = Mall_L2_norm_numerator_sum + (abs(Mall1(i, 1) - Mall2((2*i)-1, 1)))^2;
    Mall_L2_norm_denominator_sum = Mall_L2_norm_denominator_sum + (abs(Mall2((2*i)-1, 1)))^2;
    Hall_L2_norm_numerator_sum = Hall_L2_norm_numerator_sum + (abs(Hall1(i, 1) - Hall2((2*i)-1, 1)))^2;
    Hall_L2_norm_denominator_sum = Hall_L2_norm_denominator_sum + (abs(Hall2((2*i)-1, 1)))^2;    
    
    % finding the Linf norm
    new_umax_numerator = abs(Uall1(i, 1) - Uall2((2*i)-1, 1));
    new_nmax_numerator = abs(Nall1(i, 1) - Nall2((2*i)-1, 1));
    new_mmax_numerator = abs(Mall1(i, 1) - Mall2((2*i)-1, 1));
    new_hmax_numerator = abs(Hall1(i, 1) - Hall2((2*i)-1, 1));
    new_umax_denominator = abs(Uall2((2*i)-1, 1));
    new_nmax_denominator = abs(Nall2((2*i)-1, 1));
    new_mmax_denominator = abs(Mall2((2*i)-1, 1));
    new_hmax_denominator = abs(Hall2((2*i)-1, 1));
    
    % find Linf norm
    if new_umax_numerator > Uall_max_numerator
        Uall_max_numerator = new_umax_numerator;
    end
    if new_nmax_numerator > Nall_max_numerator
        Nall_max_numerator = new_nmax_numerator;
    end
    if new_mmax_numerator > Mall_max_numerator
        Mall_max_numerator = new_mmax_numerator;
    end
    if new_hmax_numerator > Hall_max_numerator
        Hall_max_numerator = new_hmax_denominator;
    end
    
    if new_umax_denominator > Uall_max_denominator
        Uall_max_denominator = new_umax_denominator;
    end
    if new_nmax_denominator > Nall_max_denominator
        Nall_max_denominator = new_nmax_denominator;
    end
    if new_mmax_denominator > Mall_max_denominator
        Mall_max_denominator = new_mmax_denominator;
    end
    if new_hmax_denominator > Hall_max_denominator
        Hall_max_denominator = new_hmax_denominator;
    end
    
end

% Linf norms
Linf_Uall = round(vpa(Uall_max_numerator/Uall_max_denominator), 5)
Linf_Nall = round(vpa(Nall_max_numerator/Nall_max_denominator), 5)
Linf_Mall = round(vpa(Mall_max_numerator/Mall_max_denominator), 5)
Linf_Hall = round(vpa(Hall_max_numerator/Hall_max_denominator), 5)

% L1 norms
L1_Uall = round(vpa(Uall_L1_norm_numerator/Uall_L1_norm_denominator), 5)
L1_Nall = round(vpa(Nall_L1_norm_numerator/Nall_L1_norm_denominator), 5)
L1_Mall = round(vpa(Mall_L1_norm_numerator/Mall_L1_norm_denominator), 5)
L1_Hall = round(vpa(Hall_L1_norm_numerator/Hall_L1_norm_denominator), 5)

% L2 norms
L2_Uall = round(vpa(sqrt(Uall_L2_norm_numerator_sum)/sqrt(Uall_L2_norm_denominator_sum)), 5)
L2_Nall = round(vpa(sqrt(Nall_L2_norm_numerator_sum)/sqrt(Nall_L2_norm_denominator_sum)), 5)
L2_Mall = round(vpa(sqrt(Mall_L2_norm_numerator_sum)/sqrt(Mall_L2_norm_denominator_sum)), 5)
L2_Hall = round(vpa(sqrt(Hall_L2_norm_numerator_sum)/sqrt(Hall_L2_norm_denominator_sum)), 5)


% compute the Linf norm difference of -all2s and -all3s
% compute the L1 norm difference of -all2s and -all3s
% compute the L2 norm difference of -all2s and -all3s
Uall_max_numerator = 0;
Nall_max_numerator = 0;
Mall_max_numerator = 0; 
Hall_max_numerator = 0;
Uall_max_denominator = 0;
Nall_max_denominator = 0;
Mall_max_denominator = 0; 
Hall_max_denominator = 0;
Uall_L1_norm_numerator = 0;
Uall_L1_norm_denominator = 0;
Uall_L2_norm_numerator_sum = 0;
Uall_L2_norm_denominator_sum = 0;
Nall_L1_norm_numerator = 0;
Nall_L1_norm_denominator = 0;
Nall_L2_norm_numerator_sum = 0;
Nall_L2_norm_denominator_sum = 0;
Mall_L1_norm_numerator = 0;
Mall_L1_norm_denominator = 0;
Mall_L2_norm_numerator_sum = 0;
Mall_L2_norm_denominator_sum = 0;
Hall_L1_norm_numerator = 0;
Hall_L1_norm_denominator = 0;
Hall_L2_norm_numerator_sum = 0;
Hall_L2_norm_denominator_sum = 0;

for i = 1:Uall1_num_of_rows
    
    % finding the L1 norm
    Uall_L1_norm_numerator = Uall_L1_norm_numerator + abs(Uall2((2*i)-1, 1) - Uall3(10*(i-1)+1, 1));
    Uall_L1_norm_denominator = Uall_L1_norm_denominator + abs(Uall3(10*(i-1)+1, 1));
    Nall_L1_norm_numerator = Nall_L1_norm_numerator + abs(Nall2((2*i)-1, 1) - Nall3(10*(i-1)+1, 1));
    Nall_L1_norm_denominator = Nall_L1_norm_denominator + abs(Nall3(10*(i-1)+1, 1));
    Mall_L1_norm_numerator = Mall_L1_norm_numerator + abs(Mall2((2*i)-1, 1) - Mall3(10*(i-1)+1, 1));
    Mall_L1_norm_denominator = Mall_L1_norm_denominator + abs(Mall3(10*(i-1)+1, 1));
    Hall_L1_norm_numerator = Hall_L1_norm_numerator + abs(Hall2((2*i)-1, 1) - Hall3(10*(i-1)+1, 1));
    Hall_L1_norm_denominator = Hall_L1_norm_denominator + abs(Hall3(10*(i-1)+1, 1));
    
    
    % finding the L2 norm
    Uall_L2_norm_numerator_sum = Uall_L2_norm_numerator_sum + (abs(Uall2((2*i)-1, 1) - Uall3(10*(i-1)+1, 1)))^2;
    Uall_L2_norm_denominator_sum = Uall_L2_norm_denominator_sum + (abs(Uall3(10*(i-1)+1, 1)))^2;
    Nall_L2_norm_numerator_sum = Nall_L2_norm_numerator_sum + (abs(Nall2((2*i)-1, 1) - Nall3(10*(i-1)+1, 1)))^2;
    Nall_L2_norm_denominator_sum = Nall_L2_norm_denominator_sum + (abs(Nall3(10*(i-1)+1, 1)))^2;
    Mall_L2_norm_numerator_sum = Mall_L2_norm_numerator_sum + (abs(Mall2((2*i)-1, 1) - Mall3(10*(i-1)+1, 1)))^2;
    Mall_L2_norm_denominator_sum = Mall_L2_norm_denominator_sum + (abs(Mall3(10*(i-1)+1, 1)))^2;
    Hall_L2_norm_numerator_sum = Hall_L2_norm_numerator_sum + (abs(Hall2((2*i)-1, 1) - Hall3(10*(i-1)+1, 1)))^2;
    Hall_L2_norm_denominator_sum = Hall_L2_norm_denominator_sum + (abs(Hall3(10*(i-1)+1, 1)))^2;    
    
    
    % finding the Linf norm
    new_umax_numerator = abs(Uall2((2*i)-1, 1) - Uall3(10*(i-1)+1, 1));
    new_nmax_numerator = abs(Nall2((2*i)-1, 1) - Nall3(10*(i-1)+1, 1));
    new_mmax_numerator = abs(Mall2((2*i)-1, 1) - Mall3(10*(i-1)+1, 1));
    new_hmax_numerator = abs(Hall2((2*i)-1, 1) - Hall3(10*(i-1)+1, 1));
    new_umax_denominator = abs(Uall3(10*(i-1)+1, 1));
    new_nmax_denominator = abs(Nall3(10*(i-1)+1, 1));
    new_mmax_denominator = abs(Mall3(10*(i-1)+1, 1));
    new_hmax_denominator = abs(Hall3(10*(i-1)+1, 1));
    
    % find Linf norm
    if new_umax_numerator > Uall_max_numerator
        Uall_max_numerator = new_umax_numerator;
    end
    if new_nmax_numerator > Nall_max_numerator
        Nall_max_numerator = new_nmax_numerator;
    end
    if new_mmax_numerator > Mall_max_numerator
        Mall_max_numerator = new_mmax_numerator;
    end
    if new_hmax_numerator > Hall_max_numerator
        Hall_max_numerator = new_hmax_denominator;
    end
    
    if new_umax_denominator > Uall_max_denominator
        Uall_max_denominator = new_umax_denominator;
    end
    if new_nmax_denominator > Nall_max_denominator
        Nall_max_denominator = new_nmax_denominator;
    end
    if new_mmax_denominator > Mall_max_denominator
        Mall_max_denominator = new_mmax_denominator;
    end
    if new_hmax_denominator > Hall_max_denominator
        Hall_max_denominator = new_hmax_denominator;
    end
    
end

% Linf norms
Linf_Uall = round(vpa(Uall_max_numerator/Uall_max_denominator), 5)
Linf_Nall = round(vpa(Nall_max_numerator/Nall_max_denominator), 5)
Linf_Mall = round(vpa(Mall_max_numerator/Mall_max_denominator), 5)
Linf_Hall = round(vpa(Hall_max_numerator/Hall_max_denominator), 5)

% L1 norms
L1_Uall = round(vpa(Uall_L1_norm_numerator/Uall_L1_norm_denominator), 5)
L1_Nall = round(vpa(Nall_L1_norm_numerator/Nall_L1_norm_denominator), 5)
L1_Mall = round(vpa(Mall_L1_norm_numerator/Mall_L1_norm_denominator), 5)
L1_Hall = round(vpa(Hall_L1_norm_numerator/Hall_L1_norm_denominator), 5)

% L2 norms
L2_Uall = round(vpa(sqrt(Uall_L2_norm_numerator_sum)/sqrt(Uall_L2_norm_denominator_sum)), 5)
L2_Nall = round(vpa(sqrt(Nall_L2_norm_numerator_sum)/sqrt(Nall_L2_norm_denominator_sum)), 5)
L2_Mall = round(vpa(sqrt(Mall_L2_norm_numerator_sum)/sqrt(Mall_L2_norm_denominator_sum)), 5)
L2_Hall = round(vpa(sqrt(Hall_L2_norm_numerator_sum)/sqrt(Hall_L2_norm_denominator_sum)), 5)








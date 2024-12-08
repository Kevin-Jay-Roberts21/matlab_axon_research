clear all
close all
clc

% Loading in Squid Axon data
data1 = load('dimnal_stim_0.003912.mat');
data2 = load('dimnal_stim_0.003911.mat');
data3 = load('dimnal_stim_0.05.mat');
data4 = load('dimless_stim_0.003912.mat');
data5 = load('dimless_stim_0.003911.mat');
data6 = load('dimless_stim_0.05.mat');

% Loading in Myelinated Axon data
data7 = load('myelin_stim_0.371.mat');
data8 = load('myelin_stim_0.5.mat');
data9 = load('myelin_stim_0.5_nodel_0.0005.mat')
data10 = load('myelin_stim_0.5_length_1..mat')
data11 = load('myelin_stim_0.5_length_1_2.mat')
data12 = load('myelin_stim_0.5_gNa_1.mat')
data13 = load('myelin_stim_0.5_radius_0.00005.mat')
data14 = load('myelin_stim_0.5_radius_0.00005_longer_NA0.5.mat')
data15 = load('myelin_stim_0.5_radius_0.00005_longer_NA0.4.mat')
data16 = load('myelin_stim_0.5_radius_0.00005_longer_NA0.3.mat')
data17 = load('myelin_stim_0.5_radius_0.00005_diff_equ.mat')
data18 = load('myelin_stim_0.5_radius_0.00005_starting_at_equil_15_node.mat')
data19 = load('myelin_stim_0.5_radius_0.00005_starting_at_equil_15_node_h_.0001.mat')
data20 = load('myelin_stim_0.5_radius_0.00005_starting_at_equil_15_node_h_.00005.mat')

% Loading SC data
data = load('First_SC_Code.mat');
Vm_all = data.Vm_all;
Vmy_all = data.Vmy_all;
Vm_minus_Vmy = data.Vm_minus_Vmy;
N_all = data.N_all;
M_all = data.M_all;
H_all = data.H_all;


%%%%%%%%%%%%%%%%%%%%
% PLOTTING 2 LINES %
%%%%%%%%%%%%%%%%%%%%

% d = data1.d;
% T = data1.T;
% m = data1.m;
% n = data1.n;
% k = data1.k;
% h = data1.h;
% t_c = data1.t_c;
% x_c = data1.x_c;
% 
% observed_time = 10; % in ms
% observed_position = 2; % in cm
%
% % SPATIAL PROFILE % 
% figure(1)
% t1 = linspace(0, d, m);
% plot(t1, data1.Uall(round(observed_time/k),:))
% hold on
% plot(t1, data2.Uall(round(observed_time/k),:))
% legend({sprintf('Voltage of the axon at time t = %g ms', observed_time), sprintf('$\\tilde{V_m}$ at $\\tilde{t} = %g$', observed_time/t_c)}, 'Interpreter','latex')
% ylabel("Voltage in millivolts.")
% xlabel("Length of the axon in cm.")
%
% % TEMPORAL PROFILE %
% figure(2)
% t2 = linspace(0, T, n);
% plot(t2, data1.Uall(:,round(observed_position/h)))
% hold on
% plot(t2, data2.Uall(:,round(observed_position/h)))
% legend({sprintf('Voltage at x = %g cm', observed_position), sprintf('$\\tilde{V_m}$ at $\\tilde{x} = %g$', observed_position/x_c)}, 'Interpreter','latex')
% ylabel("Voltage in millivolts.")
% xlabel("Time in milliseconds.")




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MYELINATED AXON ANIMATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d = data10.L;
% T = data10.T;
% m = data10.m;
% n = data10.n;
% k = data10.k;
% h = data10.h;
% % 
% % % SPATIAL PROFILE %
% % % x axis is the axon length
% t1 = linspace(0, d, m); 
% 
% figure(1);
% hold on;
% 
% xmin = 0;
% xmax = d;
% ymin = -90;
% ymax = 60;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Axon length in cm');
% ylabel('Voltage of axon in mV');
% 
% % Loop through each vector and plot them one by one
% for i = 1:n
%     x1 = data10.Uall(i,:);
% %     x2 = data11.Uall(i,:);
%     
%     % Plot the vector
%     plot(t1, x1, 'b-');
%     hold on
% %     plot(t1, x2, 'r-');
%     
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Time: %.2f ms', round(i*k, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
%     
%     % Add a pause to create animation effect
%     pause(0.01);
%     
%     cla;
% end

% TEMPORAL PROFILE %
% x axis is the axon time

% defining plot linspace interval
% T1 = 5;
% T2 = 6.2;
% steps = (T2-T1)/k+1;
% t2 = linspace(T1, T2, steps); 
% 
% figure(2);
% hold on;
% 
% xmin = T1;
% xmax = T2;
% ymin = -26;
% ymax = 10;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Time in ms');
% ylabel('Voltage of axon in mV');
% 
% % defining space interval
% P1 = 0.4; % (in cm)
% P2 = 0.6; % (in cm)
% 
% nodal_regions = data10.nodal_regions;
% 
% color = 'r-';
% % Loop through each vector and plot them one by one
% for i = P1/h:P2/h
%     x1 = data10.Uall(T1/k:T2/k,i);  
% %     x2 = data11.Uall(:,i);
%     
%     for j = 1:85
%         if (round(i*h, 4) == round(nodal_regions(1, j), 4)) || (round(i*h,4) == round(nodal_regions(2, j), 4))
%             color = 'b-';
%             break;
%         else
%             color = 'r-';
%         end
%     end
%     % Plot the vector
%     plot(t2, x1, color, 'LineWidth', 3);
%     hold on
% %     plot(t2, x2, 'r-');
% 
%     text(xmax - 0.25 * (xmax - xmin), ymin + 0.1 * (ymax - ymin), sprintf('Space: %.4f cm', round(i*h, 4)), 'FontSize', 12, 'BackgroundColor', 'w');
% 
%     % Add a pause to create animation effect
%     pause(0.05);
% 
%     cla;
% end






%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING ANIMATION %
%%%%% SQUID AXON %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% L = data20.L;
% T = data20.T;
% m = data20.m;
% n = data20.n;
% k = data20.k;
% h = data20.h;
% 
% % % SPATIAL PROFILE %
% % % x axis is the axon length
% t1 = linspace(0, L, m); 
% 
% figure(1);
% hold on;
% 
% xmin = 0;
% xmax = L;
% ymin = -90;
% ymax = 60;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Axon length in cm');
% ylabel('Voltage of axon in mV');
% 
% % Loop through each vector and plot them one by one
% for i = 1:n
%     x1 = data18.Uall(i,:);
%     x2 = data19.Uall(i,:);
%     x3 = data20.Uall(i,:);
% 
%     % Plot the vector
%     plot(t1, x1, 'b-');
%     hold on
%     plot(t1, x2, 'r-');
%     hold on
%     plot(t1, x3, 'b-');
% 
% 
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Time: %.3f ms', round(i*k, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
% 
%     % Add a pause to create animation effect
%     pause(0.001);
% 
%     cla;
% end

% TEMPORAL PROFILE %
% x axis is the axon time
% t2 = linspace(0, T, n); 
% 
% figure(2);
% hold on;
% 
% xmin = 0;
% xmax = T;
% ymin = -60;
% ymax = 90;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Time in ms');
% ylabel('Voltage of axon in mV');
% 
% % Loop through each vector and plot them one by one
% for i = 1:m
%     x1 = data20.Uall(:,i);  
%     % x2 = data19.Uall(:,i);
%     % x3 = data20.Uall(:,i);
% 
%     % Plot the vector
%     plot(t2, x1, 'b-');
%     hold on
%     % plot(t2, x2, 'r-');
%     % hold on
%     % plot(t2, x3, 'g-');
%     % hold on
% 
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.05 * (ymax - ymin), sprintf('Space: %.5f cm', round(i*h, 5)), 'FontSize', 12, 'BackgroundColor', 'w');
% 
%     % Add a pause to create animation effect
%     pause(0.01);
% 
%     cla;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING DIFFERENT RESOLUTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = data18.L;
% L1 = data19.L;
% L2 = data20.L;
% T = data18.T;
% T1 = data19.T;
% T2 = data20.T;
% m = data18.m;
% m1 = data19.m;
% m2 = data20.m;
% n = data18.n;
% n1 = data19.n;
% n2 = data20.n;
% k = data18.k;
% k1 = data19.k;
% k2 = data20.k;
% h = data18.h;
% h1 = data19.h;
% h2 = data20.h;
% 
% % % SPATIAL PROFILE %
% % % x axis is the axon length

% figure(1);
% hold on;
% 
% xmin = 0;
% xmax = L;
% ymin = -90;
% ymax = 60;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Axon length in cm');
% ylabel('Voltage of axon in mV');
% 
% x1 = data18.Uall(1,:);
% x2 = data19.Uall(1,:);
% 
% % Loop through each vector and plot them one by one
% for i = 1:n2
% 
%     x1 = data18.Uall(i,:);    
%     x2 = data19.Uall(i,:);
%     x3 = data20.Uall(i,:);
% 
%     x_range = [0, L];
%     t1 = linspace(x_range(1), x_range(2), length(x1));
%     t2 = linspace(x_range(1), x_range(2), length(x2));
%     t3 = linspace(x_range(1), x_range(2), length(x3));
% 
%     % Plot the vector
%     plot(t1, x1, 'b-');
%     hold on
%     plot(t2, x2, 'r-');
%     hold on
%     plot(t3, x3, 'g-');
%     hold on
% 
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.05 * (ymax - ymin), sprintf('Space: %.2f ms', round(i*k2, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
% 
%     xlim(x_range);
% 
%     % Add a pause to create animation effect
%     pause(0.01);
% 
%     cla;
% end

% TEMPORAL PROFILE %
% x axis is the axon time
% figure(2);
% hold on;
% 
% xmin = 0;
% xmax = T;
% ymin = -60;
% ymax = 90;
% 
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Time in ms');
% ylabel('Voltage of axon in mV');
% 
% x1 = data18.Uall(:,1);
% x2 = data19.Uall(:,1);
% 
% % Loop through each vector and plot them one by one
% for i = 1:m2
% 
%     if mod(i, 10) == 0
%         x1 = data18.Uall(:,round(i/10));
%     end
%     if mod(i, 2) == 0
%         x2 = data19.Uall(:,round(i/2));
%     end
%     x3 = data20.Uall(:,i);
% 
%     x_range = [0, T];
%     t1 = linspace(x_range(1), x_range(2), length(x1));
%     t2 = linspace(x_range(1), x_range(2), length(x2));
%     t3 = linspace(x_range(1), x_range(2), length(x3));
% 
%     % Plot the vector
%     plot(t1, x1, 'b-');
%     hold on
%     plot(t2, x2, 'r-');
%     hold on
%     plot(t3, x3, 'g-');
%     hold on
% 
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.05 * (ymax - ymin), sprintf('Space: %.5f cm', round(i*h2, 5)), 'FontSize', 12, 'BackgroundColor', 'w');
% 
%     xlim(x_range);
% 
%     % Add a pause to create animation effect
%     pause(0.01);
% 
%     cla;
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MULTIPLE LINES IN MYELINATED SECTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = data18.L;
% T = data18.T;
% m = data18.m;
% n = data18.n;
% k = data18.k;
% h = data18.h;
% % 
% % % now pick a position to plot all of the voltages (multiply by 10000 to get
% % % units in um)
% node1 = 0.3000;
% node2 = 0.3125;
% segments = (node2-node1)/h;
% recorded = 5;
% 
% 
% position1 = node1; % in cm
% position2 = node1 + 5*h; % in cm
% position3 = node1 + (5*h*2); % in cm
% position4 = node1 + (5*h*3); % in cm
% position5 = node1 + (5*h*4); % in cm
% position6 = node2; % in cm
% position7 = node2 + (5*h*1); % in cm
% position8 = node2 + (5*h*2); % in cm
% position9 = node2 + (5*h*3); % in cm
% position10 = node2 + (5*h*4); % in cm
% position11 = node2 + (5*h*5); % in cm
% position12 = node2 + (5*h*6); % in cm
% 
% list_of_positions = [position1
%                      position2
%                      position3
%                      position4
%                      position5
%                      position6
%                      position7
%                      position8
%                      position9
%                      position10
%                      position11
%                      position12];
% 
% % Times to observe the voltage along the axon
% time1 = 1; % in ms
% time2 = 1.5; % in ms
% time3 = 2; % in ms
% time4 = 4; % in ms
% time5 = 5; % in ms
% time6 = 6; % in ms
% time7 = 7; % in ms
% 
% list_of_times = [time1
%                  time2
%                  time3
%                  time4
%                  time5
%                  time6];
% 
% % % plotting Voltage vs Axon length
% figure(1)
% t1 = linspace(0, L, m);
% plot(t1, data18.Uall(round(time1/k),:))
% for i = 2:length(list_of_times)
%     hold on
%     plot(t1, data18.Uall(round(list_of_times(i)/k),:))
% end
% 
% % describing plots using legends
% legendStrings1 = {};
% for i  = 1:length(list_of_times)
%     legendStrings1{end+1} = sprintf('Voltage of the axon at time t = %g ms', list_of_times(i));
% end
% legend(legendStrings1, 'Interpreter','latex')
% ylabel("Voltage in millivolts.")
% xlabel("Length of the axon in cm.")
% 
% % plotting Voltage vs Time
% figure(2)
% 
% t2 = linspace(0, T, n); % FULL MATRIX
% % t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% % t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
% plot(t2, data14.Uall(:,round(position1/h)))
% for i = 2:length(list_of_positions)
%     hold on
%     plot(t2, data14.Uall(:,round(list_of_positions(i)/h)))
% end
% 
% % describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     message = "Voltage at x = %g cm";
%     if list_of_positions(i) == position1
%         message = "NODE 1: x = %g cm";
%     elseif list_of_positions(i) == position6
%         message = "NODE 2: x = %g cm";
%     else
%         message = sprintf('V at %d/5 of myelin', i-1);
%     end
% 
%     legendStrings2{end+1} = sprintf(message, list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex')
% ylabel("Voltage in millivolts.")
% xlabel("Time in milliseconds.")
% 
% % plotting N, M, H probability vs time (at a certain position)
% figure(3)
% plot(t2, data14.Nall(:,round(position3/h)))
% hold on
% plot(t2, data14.Mall(:,round(position3/h)))
% hold on
% plot(t2, data14.Hall(:,round(position3/h)))
% legendStrings3 = {
%     sprintf('N at x = %g cm', position3), ...
%     sprintf('M at x = %g cm', position3), ...
%     sprintf('H at x = %g cm', position3)};
% legend(legendStrings3, 'Interpreter','latex')
% ylabel("Probabilities of ion channels opening/closing.")
% xlabel("Time in milliseconds.")


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MULTIPLE LINES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = data12.L;
% T = data12.T;
% m = data12.m;
% n = data12.n;
% k = data12.k;
% h = data12.h;
% 
% position1 = 0.2; % in cm
% position2 = 0.3; % in cm
% position3 = 0.5; % in cm
% position4 = 0.6; % in cm
% position5 = 0.7; % in cm
% position6 = 0.8; % in cm
% position7 = 1; % in cm
% 
% list_of_positions = [position1
%                      position2
%                      position3
%                      position4
%                      position5
%                      position6];
% 
% % Times to observe the voltage along the axon
% time1 = 2; % in ms
% time2 = 3; % in ms
% time3 = 3.5; % in ms
% time4 = 4; % in ms
% time5 = 4.5; % in ms
% time6 = 5; % in ms
% time7 = 5.5; % in ms
% 
% list_of_times = [time1
%                  time2
%                  time3
%                  time4
%                  time5
%                  time6];
% 
% % plotting Voltage vs Axon length
% figure(1)
% t1 = linspace(0, L, m);
% plot(t1, data12.Uall(round(time1/k),:))
% for i = 2:length(list_of_times)
%     hold on
%     plot(t1, data12.Uall(round(list_of_times(i)/k),:))
% end
% 
% % describing plots using legends
% legendStrings1 = {};
% for i  = 1:length(list_of_times)
%     legendStrings1{end+1} = sprintf('Voltage of the axon at time t = %g ms', list_of_times(i));
% end
% legend(legendStrings1, 'Interpreter','latex')
% ylabel("Voltage in millivolts.")
% xlabel("Length of the axon in cm.")
% 
% % plotting Voltage vs Time
% figure(2)
% 
% t2 = linspace(0, T, n); % FULL MATRIX
% % t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% % t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
% plot(t2, data12.Uall(:,round(position1/h)))
% for i = 2:length(list_of_positions)
%     hold on
%     plot(t2, data12.Uall(:,round(list_of_positions(i)/h)))
% end
% 
% % describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     message = "Voltage at x = %g cm";
%     legendStrings2{end+1} = sprintf(message, list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex')
% ylabel("Voltage in millivolts.")
% xlabel("Time in milliseconds.")
% 
% % plotting N, M, H probability vs time (at a certain position)
% figure(3)
% plot(t2, data12.Nall(:,round(position3/h)))
% hold on
% plot(t2, data12.Mall(:,round(position3/h)))
% hold on
% plot(t2, data12.Hall(:,round(position3/h)))
% legendStrings3 = {
%     sprintf('N at x = %g cm', position3), ...
%     sprintf('M at x = %g cm', position3), ...
%     sprintf('H at x = %g cm', position3)};
% legend(legendStrings3, 'Interpreter','latex')
% ylabel("Probabilities of ion channels opening/closing.")
% xlabel("Time in milliseconds.")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING SC CODE RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = data.L;
% T = data.T;
% n = data.n;
% m = data.m;
% dx = data.dx;
% dt = data.dt;
% 
% position1 = 0.01; % in cm
% position2 = 0.015; % in cm
% position3 = 0.02; % in cm
% position4 = 0.025; % in cm
% position5 = 0.03; % in cm
% position6 = 0.035; % in cm
% position7 = 0.04; % in cm
% 
% list_of_positions = [position1
%                      position7];
% 
% % Times to observe the voltage along the axon
% time1 = 1; % in ms
% time2 = 2; % in ms
% time3 = 4; % in ms
% time4 = 5; % in ms
% time5 = 7; % in ms
% time6 = 8; % in ms
% time7 = 9; % in ms
% 
% list_of_times = [time1
%                  time2
%                  time3
%                  time4
%                  ];
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOTTING Vm SPATIAL AND TEMPORAL PROFILES %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First figure: Voltage along the axon at different times
% figure(1);
% 
% % Adjust the figure size (Position [left, bottom, width, height])
% set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)
% 
% % Create subplot (1 row, 2 columns, 1st subplot)
% subplot(1, 2, 1);
% t1 = linspace(0, 10000*L, m);
% plot(t1, Vm_all(round(time1/dt),:));
% hold on
% for i = 2:length(list_of_times)
%     plot(t1, Vm_all(round(list_of_times(i)/dt),:));
% end
% 
% % Describing plots using legends
% legendStrings1 = {};
% for i  = 1:length(list_of_times)
%     legendStrings1{end+1} = sprintf('$V_m$ at t = %g ms', list_of_times(i));
% end
% legend(legendStrings1, 'Interpreter','latex');
% ylabel("$V_m$ in millivolts.", 'Interpreter','latex');
% xlabel("Length of the axon in um.");
% 
% 
% % Second figure: Voltage vs Time at different positions
% % Create subplot (1 row, 2 columns, 2nd subplot)
% subplot(1, 2, 2);
% t2 = linspace(0, T, n); % FULL MATRIX
% plot(t2, Vm_all(:,round(position1/dx)));
% hold on
% for i = 2:length(list_of_positions)
%     plot(t2, Vm_all(:,round(list_of_positions(i)/dx)));
% end
% 
% % Describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     legendStrings2{end+1} = sprintf('$V_m$ at x = %g um', 10000*list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex');
% ylabel("$V_m$ in millivolts.", 'Interpreter', 'latex');
% xlabel("Time in milliseconds.");
% 
% 
% % PLOTTING Vmy SPATIAL AND TEMPORAL PROFILES 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First figure: Voltage along the axon at different times
% figure(2);
% 
% % Adjust the figure size (Position [left, bottom, width, height])
% set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)
% 
% % Create subplot (1 row, 2 columns, 1st subplot)
% subplot(1, 2, 1);
% t1 = linspace(0, L*10000, m);
% plot(t1, Vmy_all(round(time1/dt),:));
% hold on
% for i = 2:length(list_of_times)
%     plot(t1, Vmy_all(round(list_of_times(i)/dt),:));
% end
% 
% % Describing plots using legends
% legendStrings1 = {};
% for i  = 1:length(list_of_times)
%     legendStrings1{end+1} = sprintf('$V_{my}$ at t = %g ms', list_of_times(i));
% end
% legend(legendStrings1, 'Interpreter','latex');
% ylabel("$V_{my}$ in millivolts.", 'Interpreter','latex');
% xlabel("Length of the axon in um.");
% 
% 
% % Second figure: Voltage vs Time at different positions
% % Create subplot (1 row, 2 columns, 2nd subplot)
% subplot(1, 2, 2);
% t2 = linspace(0, T, n); % FULL MATRIX
% plot(t2, Vmy_all(:,round(position1/dx)));
% hold on
% for i = 2:length(list_of_positions)
%     plot(t2, Vmy_all(:,round(list_of_positions(i)/dx)));
% end
% 
% % Describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     legendStrings2{end+1} = sprintf('$V_{my}$ at x = %g um', 10000*list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex');
% ylabel("$V_{my}$ in millivolts.", 'Interpreter', 'latex');
% xlabel("Time in milliseconds.");
% 
% 
% % PLOTTING Vm-Vmy SPATIAL AND TEMPORAL PROFILES 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First figure: Voltage along the axon at different times
% figure(3);
% 
% % Adjust the figure size (Position [left, bottom, width, height])
% set(gcf, 'Position', [100, 100, 1200, 500]); % Increase width (1200)
% 
% % Create subplot (1 row, 2 columns, 1st subplot)
% subplot(1, 2, 1);
% t1 = linspace(0, 10000*L, m);
% plot(t1, Vm_minus_Vmy(round(time1/dt),:));
% hold on
% for i = 2:length(list_of_times)
%     plot(t1, Vm_minus_Vmy(round(list_of_times(i)/dt),:));
% end
% 
% % Describing plots using legends
% legendStrings1 = {};
% for i  = 1:length(list_of_times)
%     legendStrings1{end+1} = sprintf('$V_m - V_{my}$ at t = %g ms', list_of_times(i));
% end
% legend(legendStrings1, 'Interpreter','latex');
% ylabel("$V_m - V_{my}$ in millivolts.", 'Interpreter','latex');
% xlabel("Length of the axon in um.");
% 
% 
% % Second figure: Voltage vs Time at different positions
% % Create subplot (1 row, 2 columns, 2nd subplot)
% subplot(1, 2, 2);
% t2 = linspace(0, T, n); % FULL MATRIX
% plot(t2, Vm_minus_Vmy(:,round(position1/dx)));
% hold on
% for i = 2:length(list_of_positions)
%     plot(t2, Vm_minus_Vmy(:,round(list_of_positions(i)/dx)));
% end
% 
% % Describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     legendStrings2{end+1} = sprintf('$V_m - V_{my}$ at x = %g um', 10000*list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex');
% ylabel("$V_m - V_{my}$ in millivolts.", 'Interpreter', 'latex');
% xlabel("Time in milliseconds.");
% 
% 
% % PLOTTING n, m, and h TEMPROAL PROFILES 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% plot(t2, N_all(:,round(position3/dx)))
% hold on
% plot(t2, M_all(:,round(position3/dx)))
% hold on
% plot(t2, H_all(:,round(position3/dx)))
% legendStrings3 = {
%     sprintf('n at x = %g um', 10000*position3), ...
%     sprintf('m at x = %g um', 10000*position3), ...
%     sprintf('h at x = %g um', 10000*position3)};
% legend(legendStrings3, 'Interpreter','latex')
% ylabel("Probabilities of ion channels opening/closing.")
% xlabel("Time in milliseconds.")

clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%
% PLOTTING 2 LINES %
%%%%%%%%%%%%%%%%%%%%

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


% norms of each matrix comparison
% conclude that the matrices are very similar due to 10^-9 norm or smaller
norm(data1.Uall-data4.Uall*25)
norm(data2.Uall-data5.Uall*25)
norm(data3.Uall-data6.Uall*25)

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MYELINATED AXON %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d = data8.d;
% T = data8.T;
% m = data8.m;
% n = data8.n;
% k = data8.k;
% h = data8.h;

% % SPATIAL PROFILE %
% % x axis is the axon length
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
%     x1 = data8.Uall(i,:);
% %     x2 = data2.Uall(i,:);
%     
%     % Plot the vector
%     plot(t1, x1, 'b-');
% %     hold on
% %     plot(t1, x2, 'r-');
%     
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Time: %.2f ms', round(i*k, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
%     
%     % Add a pause to create animation effect
%     pause(0.01);
%     
%     cla;
% end

% % TEMPORAL PROFILE %
% x axis is the axon time
% t2 = linspace(0, T, n); 
% 
% figure(2);
% hold on;
% 
% xmin = 0;
% xmax = T;
% ymin = -90;
% ymax = 60;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Time in ms');
% ylabel('Voltage of axon in mV');
% 
% % Loop through each vector and plot them one by one
% for i = 1:m
%     x1 = data8.Uall(:,i);  
% %     x2 = data8.Uall(:,i);
%     
%     % Plot the vector
%     plot(t2, x1, 'b-');
%     hold on
% %     plot(t2, x2, 'r-');
%     
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Space: %.2f cm', round(i*h, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
%     
%     % Add a pause to create animation effect
%     pause(0.01);
%     
%     cla;
% end






%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING ANIMATION %
%%%%% SQUID AXON %%%%%
%%%%%%%%%%%%%%%%%%%%%%

d = data1.d;
T = data1.T;
m = data1.m;
n = data1.n;
k = data1.k;
h = data1.h;

% SPATIAL PROFILE %
% x axis is the axon length
t1 = linspace(0, d, m); 

figure(1);
hold on;

xmin = 0;
xmax = d;
ymin = -90;
ymax = 60;

axis([xmin xmax ymin ymax]);  % Set axis limits
xlabel('Axon length in cm');
ylabel('Voltage of axon in mV');

% Loop through each vector and plot them one by one
for i = 1:n
    x1 = data1.Uall(i,:);
    x2 = data4.Uall(i,:);
    
    % Plot the vector
    plot(t1, x1, 'b-');
    hold on
    plot(t1, x2, 'r-');
    
    text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Time: %.2f ms', round(i*k, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
    
    % Add a pause to create animation effect
    pause(0.001);
    
    cla;
end
% 
% % TEMPORAL PROFILE %
% % x axis is the axon time
% t2 = linspace(0, T, n); 
% 
% figure(2);
% hold on;
% 
% xmin = 0;
% xmax = T;
% ymin = -90;
% ymax = 60;
% 
% axis([xmin xmax ymin ymax]);  % Set axis limits
% xlabel('Time in ms');
% ylabel('Voltage of axon in mV');
% 
% % Loop through each vector and plot them one by one
% for i = 1:m
%     x1 = data1.Uall(:,i);  
%     x2 = data2.Uall(:,i);
%     
%     % Plot the vector
%     plot(t2, x1, 'b-');
%     hold on
%     plot(t2, x2, 'r-');
%     
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Space: %.2f cm', round(i*h, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
%     
%     % Add a pause to create animation effect
%     pause(0.01);
%     
%     cla;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MULTIPLE LINES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % now pick a position to plot all of the voltages (multiply by 10000 to get
% % units in um)
% position1 = 0.5; % in cm
% position2 = 1; % in cm
% position3 = 1.5; % in cm
% position4 = 2; % in cm
% position5 = 2.5; % in cm
% position6 = 3; % in cm
% position7 = 3.5; % in cm
% 
% list_of_positions = [position1
%                      position2
%                      position3
%                      position4
%                      position5
%                      position6
%                      position7];
% 
% % Times to observe the voltage along the axon
% time1 = 5; % in ms
% time2 = 6; % in ms
% time3 = 8; % in ms
% time4 = 10; % in ms
% time5 = 10.5; % in ms
% time6 = 10.8; % in ms
% time7 = 11; % in ms
% 
% list_of_times = [time1
%                  time2
%                  time3
%                  time4
%                  time5
%                  time6
%                  time7];
% 
% % plotting Voltage vs Axon length
% figure(1)
% t1 = linspace(0, d, m);
% plot(t1, Uall(round(time1/k),:))
% for i = 2:length(list_of_times)
%     hold on
%     plot(t1, Uall(round(list_of_times(i)/k),:))
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
% t2 = linspace(0, T, n); % FULL MATRIX
% % t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% % t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
% plot(t2, Uall(:,round(position1/h)))
% for i = 2:length(list_of_positions)
%     hold on
%     plot(t2, Uall(:,round(list_of_positions(i)/h)))
% end
% 
% % describing plots using legends
% legendStrings2 = {};
% for i = 1:length(list_of_positions)
%     legendStrings2{end+1} = sprintf('Voltage at x = %g cm', list_of_positions(i));
% end
% legend(legendStrings2, 'Interpreter', 'latex')
% ylabel("Voltage in millivolts.")
% xlabel("Time in milliseconds.")
% 
% % plotting N, M, H probability vs time (at a certain position)
% figure(3)
% plot(t2, Nall(:,round(position3/h)))
% hold on
% plot(t2, Mall(:,round(position3/h)))
% hold on
% plot(t2, Hall(:,round(position3/h)))
% legendStrings3 = {
%     sprintf('N at x = %g cm', position3), ...
%     sprintf('M at x = %g cm', position3), ...
%     sprintf('H at x = %g cm', position3)};
% legend(legendStrings3, 'Interpreter','latex')
% ylabel("Probabilities of ion channels opening/closing.")
% xlabel("Time in milliseconds.")
% 
% 
% % save U,N,,'cable_0.01'  % this is for...
% % save(U, M, N, H, 'test' % this )
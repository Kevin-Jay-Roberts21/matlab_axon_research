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
data9 = load('myelin_stim_0.5_nodel_0.0005.mat')
data10 = load('myelin_stim_0.5_length_1..mat')
data11 = load('myelin_stim_0.5_length_1_2.mat')


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

d = data10.L;
T = data10.T;
m = data10.m;
n = data10.n;
k = data10.k;
h = data10.h;
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

% d = data1.d;
% T = data1.T;
% m = data1.m;
% n = data1.n;
% k = data1.k;
% h = data1.h;

% SPATIAL PROFILE %
% x axis is the axon length
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
%     x1 = data1.Uall(i,:);
%     % x2 = data4.Uall(i,:);
%     
%     % Plot the vector
%     plot(t1, x1, 'b-');
%     hold on
%     % plot(t1, x2, 'r-');
%     
%     text(xmin + 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin), sprintf('Time: %.2f ms', round(i*k, 2)), 'FontSize', 12, 'BackgroundColor', 'w');
%     
%     % Add a pause to create animation effect
%     pause(0.001);
%     
%     cla;
% end
% 
% TEMPORAL PROFILE %
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
%     x1 = data1.Uall(:,i);  
% %     x2 = data2.Uall(:,i);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MULTIPLE LINES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
node1 = 0.4440;
node2 = 0.4565;
segments = (node2-node1)/h;
recorded = 5;


position1 = node1; % in cm
position2 = node1 + 5*h; % in cm
position3 = node1 + (5*h*2); % in cm
position4 = node1 + (5*h*3); % in cm
position5 = node2; % in cm
position6 = 0.4560; % in cm
position7 = 0.4565; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5];

% Times to observe the voltage along the axon
time1 = 2; % in ms
time2 = 4; % in ms
time3 = 5; % in ms
time4 = 6; % in ms
time5 = 7; % in ms
time6 = 8; % in ms
time7 = 9; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];

% % plotting Voltage vs Axon length
% figure(1)
% t1 = linspace(0, d, m);
% plot(t1, data10.Uall(round(time1/k),:))
% for i = 2:length(list_of_times)
%     hold on
%     plot(t1, data10.Uall(round(list_of_times(i)/k),:))
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

% plotting Voltage vs Time
figure(2)

t2 = linspace(0, T, n); % FULL MATRIX
% t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, data10.Uall(:,round(position1/h)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, data10.Uall(:,round(list_of_positions(i)/h)))
end

% describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    message = "Voltage at x = %g cm";
    if list_of_positions(i) == 0.4440
        message = "NODE 1: x = %g cm";
    elseif list_of_positions(i) == 0.4565
        message = "NODE 2: x = %g cm";
    end
    
    legendStrings2{end+1} = sprintf(message, list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

% plotting N, M, H probability vs time (at a certain position)
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

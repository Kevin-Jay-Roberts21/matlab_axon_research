

% Comparing the dimensionless vs dimensional simulations
dimnal_stim_003948 = load('dimnal_stim_0.003948.mat');
dimnal_stim_003947 = load('dimnal_stim_0.003947.mat');
dimnal_stim_05 = load('dimnal_stim_0.05.mat');
dimless_stim_003948 = load('dimless_stim_0.003948.mat');
dimless_stim_003947 = load('dimless_stim_0.003947.mat');
dimless_stim_05 = load('dimless_stim_0.05.mat');

% norms of each matrix comparison
% conclude that the matrices are very similar due to 10^-9 norm or smaller
norm(dimnal_stim_003948.Uall-dimless_stim_003948.Uall*25)
norm(dimnal_stim_003947.Uall-dimless_stim_003947.Uall*25)
norm(dimnal_stim_05.Uall-dimless_stim_05.Uall*25)

d = dimnal_stim_003948.d;
T = dimnal_stim_003948.T;
m = dimnal_stim_003948.m;
n = dimnal_stim_003948.n;
k = dimnal_stim_003948.k;
t_c = dimless_stim_003948.t_c;
x_c = dimless_stim_003948.x_c;

observed_time = 10; % in ms
observed_position = 2; % in cm

figure(1)
t1 = linspace(0, d, m);
plot(t1, dimnal_stim_003948.Uall(round(observed_time/k),:))
hold on
plot(t1, dimless_stim_003948.Uall(round(observed_time/k),:))
legend({sprintf('Voltage of the axon at time t = %g ms', observed_time), sprintf('$\\tilde{V_m}$ at $\\tilde{t} = %g$', observed_time/t_c)}, 'Interpreter','latex')
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in cm.")

figure(2)
t2 = linspace(0, T, n);
plot(t2, dimnal_stim_003948.Uall(:,round(observed_position/h)))
hold on
plot(t2, dimless_stim_003948.Uall(:,round(observed_position/h)))
legend({sprintf('Voltage at x = %g cm', observed_position), sprintf('$\\tilde{V_m}$ at $\\tilde{x} = %g$', observed_position/x_c)}, 'Interpreter','latex')
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

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
% defining matrices from saved data
Uall = load('saved_Uall').Uall;
Nall = load('saved_Nall').Nall;
Mall = load('saved_Mall').Mall;
Hall = load('saved_Hall').Hall;


% now pick a position to plot all of the voltages (multiply by 10000 to get
% units in um)
position1 = 0.5; % in cm
position2 = 1; % in cm
position3 = 1.5; % in cm
position4 = 2; % in cm
position5 = 2.5; % in cm
position6 = 3; % in cm
position7 = 3.5; % in cm

list_of_positions = [position1
                     position2
                     position3
                     position4
                     position5
                     position6
                     position7];

% Times to observe the voltage along the axon
time1 = 5; % in ms
time2 = 6; % in ms
time3 = 8; % in ms
time4 = 10; % in ms
time5 = 10.5; % in ms
time6 = 10.8; % in ms
time7 = 11; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];

% plotting Voltage vs Axon length
figure(1)
t1 = linspace(0, d, m);
plot(t1, Uall(round(time1/k),:))
for i = 2:length(list_of_times)
    hold on
    plot(t1, Uall(round(list_of_times(i)/k),:))
end

% describing plots using legends
legendStrings1 = {};
for i  = 1:length(list_of_times)
    legendStrings1{end+1} = sprintf('Voltage of the axon at time t = %g ms', list_of_times(i));
end
legend(legendStrings1, 'Interpreter','latex')
ylabel("Voltage in millivolts.")
xlabel("Length of the axon in cm.")

% plotting Voltage vs Time
figure(2)
t2 = linspace(0, T, n); % FULL MATRIX
% t2 = linspace(0, T, n*k*2); % MATRIX AT EVERY 50th iteration
% t2 = linspace(0, T, n*k); % MATRIX AT EVERY 100th iteration
plot(t2, Uall(:,round(position1/h)))
for i = 2:length(list_of_positions)
    hold on
    plot(t2, Uall(:,round(list_of_positions(i)/h)))
end

% describing plots using legends
legendStrings2 = {};
for i = 1:length(list_of_positions)
    legendStrings2{end+1} = sprintf('Voltage at x = %g cm', list_of_positions(i));
end
legend(legendStrings2, 'Interpreter', 'latex')
ylabel("Voltage in millivolts.")
xlabel("Time in milliseconds.")

% plotting N, M, H probability vs time (at a certain position)
figure(3)
plot(t2, Nall(:,round(position3/h)))
hold on
plot(t2, Mall(:,round(position3/h)))
hold on
plot(t2, Hall(:,round(position3/h)))
legendStrings3 = {
    sprintf('N at x = %g cm', position3), ...
    sprintf('M at x = %g cm', position3), ...
    sprintf('H at x = %g cm', position3)};
legend(legendStrings3, 'Interpreter','latex')
ylabel("Probabilities of ion channels opening/closing.")
xlabel("Time in milliseconds.")


% save U,N,,'cable_0.01'  % this is for...
% save(U, M, N, H, 'test' % this )
% Axon voltage plotter. This code is meant to plot data from the following
% files: 
% squid_giant_axon.m 
% myelianted_HH_axon.m,
% dimensionless_squid_giant_axon.m 
% sc_code_implicit_v1.m
% sc_code_implicit_v2.m
% sc_code_implicit_v3.m
% dc_code_implicit.m
% 
% Additionally, variables like list_of_positions, list_of_times, and pause 
% may have to be modified depending on HH, SC or DC to obtain useful animations
% and/or plots.
% 
% Kevin Roberts
% January 2025

% Plotting HH or SC or DC code can be

clear all
close all
clc

% getting the saved data
% HH_data = load('HH_data.mat');
% HH_data1 = load('HH_data1.mat');
% HH_data2 = load('HH_data2.mat');
% SC_data = load('SC_data.mat');
% SC_data1 = load('SC_data1.mat');
% SC_data2 = load('SC_data2.mat');

HH_data_Temp_default = load('axon_simulations/HH_temp_data/HH_data_Temp_6.3.mat');
HH_data_Temp_7 = load('axon_simulations/HH_temp_data/HH_data_Temp_7.mat');
HH_data_Temp_8 = load('axon_simulations/HH_temp_data/HH_data_Temp_8.mat');
HH_data_Temp_9 = load('axon_simulations/HH_temp_data/HH_data_Temp_9.mat');
HH_data_Temp_10 = load('axon_simulations/HH_temp_data/HH_data_Temp_10.mat');
HH_data_Temp_15 = load('axon_simulations/HH_temp_data/HH_data_Temp_15.mat');
HH_data_Temp_20 = load('axon_simulations/HH_temp_data/HH_data_Temp_20.mat');
HH_data_Temp_25 = load('axon_simulations/HH_temp_data/HH_data_Temp_25.mat');
HH_data_Temp_30 = load('axon_simulations/HH_temp_data/HH_data_Temp_30.mat');
HH_data_Temp_35 = load('axon_simulations/HH_temp_data/HH_data_Temp_35.mat');

% picking time shots
time1 = 5; % in ms
time2 = 6; % in ms
time3 = 8; % in ms
time4 = 10; % in ms
time5 = 11; % in ms
time6 = 12; % in ms
time7 = 13; % in ms

list_of_times = [time1
                 time2
                 time3
                 time4
                 time5
                 time6
                 time7];


% picking space shots
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


% picking pause (this controls the speed of the animation, the pause variable 
% is in seconds. It is how many seconds each time (or space) shot will be
% paused at). NOTE: the legends are what slows down the animations, may
% condsider getting rid of or adding them in certain cases.
p = 0.001;

% creating a set of data from multiple experiments used to plot animation
% (the first element in the set_of_data is darkred, then the proceeding elements get
% brighter and brighter until the last element which is the brightest red)
set_of_data = [HH_data_Temp_default, HH_data_Temp_7, HH_data_Temp_8, HH_data_Temp_9, HH_data_Temp_10, HH_data_Temp_15, HH_data_Temp_20, HH_data_Temp_25, HH_data_Temp_30];


% plot_animation_voltage_vs_time(SC_data2, p);
% plot_animation_voltage_vs_space(HH_data_Temp_30, p);
% plot_animation_probabilities(SC_data, p);
% plot_time_and_space_shots(HH_data, list_of_positions, list_of_times);
% plot_voltage_vs_time_comparison(set_of_data, p);
plot_voltage_vs_space_comparison(set_of_data, p);




%%%%%%%%%%%%%%%%%%%%%
% PLOTTER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%

% Animation that plots the voltage vs space
function plot_animation_voltage_vs_space(data, p)
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % SPATIAL PROFILE %
    % x axis is the axon length
    t = linspace(0, L, m); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Axon length in cm');
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')

    % Loop through each vector and plot them one by one
    for i = 1:n
        x = data.Vm_all(i,:);
        
        % Plot the vector
        plot(t, x, 'b-');
        
        text(xmin + 0.05, ymax + 0.8, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

% Animation that plots the voltage vs time
function plot_animation_voltage_vs_time(data, p)
    T = data.T;
    m = data.m;
    n = data.n;
    dx = data.dx;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 90;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Time in milliseconds');
    ylabel('$V_m$ in millivolts', 'Interpreter','latex');

    % Loop through each vector and plot them one by one
    for i = 1:m
        x = data.Vm_all(:,i);  

        % Plot the vector
        plot(t, x, 'b-');
        hold on

        text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

% Animation that plots the probabilities vs time
function plot_animation_probabilities(data, p)
    T = data.T;
    m = data.m;
    n = data.n;
    dx = data.dx;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = 0;
    ymax = 1;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Time in milliseconds');
    ylabel('Probabilities of gates begin open.');
    
    % Loop through each vector and plot them one by one
    for i = 1:m
        x1 = data.N_all(:,i);  
        x2 = data.M_all(:,i);  
        x3 = data.H_all(:,i);  

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        p3 = plot(t, x3, 'g-');
        hold on
        
        % Add the legend (NOTE: legends are what slows down animation)
        legend([p1, p2, p3], 'K activation gate', 'Na activation gate', 'Na inactivation date', 'Location', 'northeast');

        text(xmin, ymax + 0.05, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');
        
        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

% plotting multiple time and space shots of voltage for squid
function plot_time_and_space_shots(data, list_of_positions, list_of_times)

    L = data.L;
    T = data.T;
    m = data.m;
    n = data.n;
    dt = data.dt;
    dx = data.dx;
    
    % plotting Voltage vs Axon length
    figure(1)
    t1 = linspace(0, L, m);
    for i = 1:length(list_of_times)
        plot(t1, data.Vm_all(round(list_of_times(i)/dt),:))
        hold on
    end
    
    % describing plots using legends
    legendStrings1 = {};
    for i  = 1:length(list_of_times)
        legendStrings1{end+1} = sprintf('$V_m$ at t = %g ms', list_of_times(i));
    end
    legend(legendStrings1, 'Interpreter','latex')
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel("Length of the axon in cm.")
    
    % plotting Voltage vs Time
    figure(2)
    t2 = linspace(0, T, n); % FULL MATRIX
    for i = 1:length(list_of_positions)
        plot(t2, data.Vm_all(:,round(list_of_positions(i)/dx)))
        hold on
    end
    
    % describing plots using legends
    legendStrings2 = {};
    for i = 1:length(list_of_positions)
        legendStrings2{end+1} = sprintf('$V_m$ at x = %g cm', list_of_positions(i));
    end
    legend(legendStrings2, 'Interpreter', 'latex')
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel("Time in milliseconds.")
    
    % plotting N, M, H probability vs time (at the first position list_of_positions(1))
    figure(3)
    plot(t2, data.N_all(:,round(list_of_positions(1)/dx)))
    hold on
    plot(t2, data.M_all(:,round(list_of_positions(1)/dx)))
    hold on
    plot(t2, data.H_all(:,round(list_of_positions(1)/dx)))
    legendStrings3 = {
        sprintf('n at x = %g cm', list_of_positions(1)), ...
        sprintf('m at x = %g cm', list_of_positions(1)), ...
        sprintf('h at x = %g cm', list_of_positions(1))};
    legend(legendStrings3, 'Interpreter','latex')
    ylabel("Probabilities of ion channels opening/closing.")
    xlabel("Time in milliseconds.")
end

% Plots animation of axon voltage in two different experiments
function plot_voltage_vs_time_comparison(data_set, p)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    T = data_set(1).T;
    m = data_set(1).m;
    n = data_set(1).n;
    dx = data_set(1).dx;
    
    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n);
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 60;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Time in milliseconds');
    ylabel('$V_m$ in millivolts', 'Interpreter','latex');

    % Loop through each vector and plot them one by one
    for i = 1:m
        
        for j = 1:length(data_set)
            plot(t, data_set(j).Vm_all(:,i), 'Color', [1, 0, 0] * j/length(data_set));
            hold on
        end
        
        text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

function plot_voltage_vs_space_comparison(data_set, p)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    L = data_set(1).L;
    m = data_set(1).m;
    n = data_set(1).n;
    dt = data_set(1).dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 60;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Length of axon in cm');
    ylabel('$V_m$ in millivolts', 'Interpreter','latex');
    
    % Loop through each vector and plot them one by one
    for i = 1:n
        for j = 1:length(data_set)
            plot(t, data_set(j).Vm_all(i,:), 'Color', [1, 0, 0] * j/length(data_set));
            hold on
        end
        
        % Add the legend (NOTE: the legend is what is slowing down the animation)
        % legend([p1, p2], 'data1 voltage', 'data2 voltage', 'Location', 'northeast');
        text(xmin + 0.2, ymax + 0.1, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end




% Axon voltage plotter. This code is meant to plot data from the following
% files: squid_giant_axon.m, myelianted_HH_axon.m,
% dimensionless_squid_giant_axon.m, and sc_code_implicit.m
% Kevin Roberts
% January 2025

% Plotting will be done in 2 main categories: HH, and SC/DC (or rat). Each 
% category will have the option to plot time and space shots, animation or
% comparisons.

clear all
close all
clc

% getting the saved data
HH_data = load('HH_data.mat');
HH_data1 = load('HH_data1.mat');
SC_data = load('First_SC_Code.mat');
% DC_data = load('dc_data.mat');

% plot_squid_animation_temporal(HH_data);
% plot_squid_animation_spatial(HH_data);
plot_squid_animation_probabilities(HH_data);
% plot_squid_time_and_space_shots(HH_data);
plot_squid_voltage_vs_time_comparison(HH_data, HH_data1);
plot_squid_voltage_vs_space_comparison(HH_data, HH_data1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SQUID PLOTTER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Animation that plots the voltage vs time
function plot_squid_animation_temporal(data)
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % % SPATIAL PROFILE %
    % % x axis is the axon length
    t = linspace(0, L, m); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 60;

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
        pause(0.001);

        cla;
    end
end

% Animation that plots the voltage vs space
function plot_squid_animation_spatial(data)
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
    ymax = 60;

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
        pause(0.01);

        cla;
    end
end

% Animation that plots the probabilities vs time
function plot_squid_animation_probabilities(data)
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
        
        % Add the legend
        legend([p1, p2, p3], 'K activation gate', 'Na activation gate', 'Na inactivation date', 'Location', 'northeast');

        text(xmin, ymax + 0.05, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');
        
        % Add a pause to create animation effect
        pause(0.01);

        cla;
    end
end

% plotting multiple time and space shots of voltage for squid
function plot_squid_time_and_space_shots(data)

    L = data.L;
    T = data.T;
    m = data.m;
    n = data.n;
    dt = data.dt;
    dx = data.dx;

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
    
    % plotting Voltage vs Axon length
    figure(1)
    t1 = linspace(0, L, m);
    plot(t1, data.Vm_all(round(time1/dt),:))
    for i = 2:length(list_of_times)
        hold on
        plot(t1, data.Vm_all(round(list_of_times(i)/dt),:))
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
    plot(t2, data.Vm_all(:,round(position1/dx)))
    for i = 2:length(list_of_positions)
        hold on
        plot(t2, data.Vm_all(:,round(list_of_positions(i)/dx)))
    end
    
    % describing plots using legends
    legendStrings2 = {};
    for i = 1:length(list_of_positions)
        legendStrings2{end+1} = sprintf('$V_m$ at x = %g cm', list_of_positions(i));
    end
    legend(legendStrings2, 'Interpreter', 'latex')
    ylabel('$V_m$ in millivolts.', 'Interpreter', 'latex')
    xlabel("Time in milliseconds.")
    
    % plotting N, M, H probability vs time (at a certain position)
    figure(3)
    plot(t2, data.N_all(:,round(position3/dx)))
    hold on
    plot(t2, data.M_all(:,round(position3/dx)))
    hold on
    plot(t2, data.H_all(:,round(position3/dx)))
    legendStrings3 = {
        sprintf('n at x = %g cm', position3), ...
        sprintf('m at x = %g cm', position3), ...
        sprintf('h at x = %g cm', position3)};
    legend(legendStrings3, 'Interpreter','latex')
    ylabel("Probabilities of ion channels opening/closing.")
    xlabel("Time in milliseconds.")
end

% Plots animation of axon voltage in two different experiments
function plot_squid_voltage_vs_time_comparison(data1, data2)

    T = data1.T;
    m = data1.m;
    n = data1.n;
    dx = data1.dx;

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
        x1 = data1.Vm_all(:,i);  
        x2 = data2.Vm_all(:,i);

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        
        % Add the legend
        legend([p1, p2], 'data1 voltage', 'data2 voltage', 'Location', 'northeast');
        text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(0.01);

        cla;
    end
end

function plot_squid_voltage_vs_space_comparison(data1, data2)
    
    L = data1.L;
    m = data1.m;
    n = data1.n;
    dt = data1.dt;

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
        x1 = data1.Vm_all(i,:);  
        x2 = data2.Vm_all(i,:);

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        
        % Add the legend
        legend([p1, p2], 'data1 voltage', 'data2 voltage', 'Location', 'northeast');
        text(xmin + 0.2, ymax + 0.1, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

        % Add a pause to create animation effect
        pause(0.01);

        cla;
    end
end



function plot_rat_animation(SC_data)
end
function plot_rat_time_and_space_shots(SC_data)
end
function plot_rat_comparison(SC_data, SC_data1)
end
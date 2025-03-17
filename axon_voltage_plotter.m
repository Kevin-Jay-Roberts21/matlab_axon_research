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
% HH_data1 = load('HH_1.mat');
% HH_data2 = load('HH_data2.mat');
% SC_data = load('SC_data.mat');
% DC_data = load('DC_data.mat');

% SC_data1 = load('projects/SC_data1.mat');
% SC_data2 = load('projects/SC_data2.mat');
% SC_data3 = load('projects/SC_data3.mat');
SC_data_Sv_200 = load('SC_data_Sv_200.mat');
DC_data_Sv_200 = load('DC_data_Sv_200.mat');


HH_data_Temp_default = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_6.3.mat');
HH_data_Temp_7 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_7.mat');
HH_data_Temp_8 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_8.mat');
HH_data_Temp_9 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_9.mat');
HH_data_Temp_10 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_10.mat');
HH_data_Temp_15 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_15.mat');
HH_data_Temp_20 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_20.mat');
HH_data_Temp_25 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_25.mat');
HH_data_Temp_30 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_30.mat');
HH_data_Temp_31 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_31.mat');
HH_data_Temp_32 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_32.mat');
HH_data_Temp_33 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_33.mat');
HH_data_Temp_34 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_34.mat');
HH_data_Temp_35 = load('projects/axon_simulations/HH_temp_data/HH_data_Temp_35.mat');

SC_data_Temp_20 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_20.mat');
SC_data_Temp_22 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_22.mat');
SC_data_Temp_24 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_24.mat');
SC_data_Temp_26 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_26.mat');
SC_data_Temp_28 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_28.mat');
SC_data_Temp_30 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_30.mat');
SC_data_Temp_32 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_32.mat');
SC_data_Temp_34 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_34.mat');
SC_data_Temp_36 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_36.mat');
SC_data_Temp_38 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_38.mat');
SC_data_Temp_40 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_40.mat');
SC_data_Temp_42 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_42.mat');
SC_data_Temp_44 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_44.mat');
SC_data_Temp_46 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_46.mat');
SC_data_Temp_48 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_48.mat');
SC_data_Temp_50 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_50.mat');
SC_data_Temp_52 = load('projects/axon_simulations/SC_temp_data/SC_data_Temp_52.mat');

DC_data_Temp_20 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_20.mat');
DC_data_Temp_22 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_22.mat');
DC_data_Temp_24 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_24.mat');
DC_data_Temp_26 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_26.mat');
DC_data_Temp_28 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_28.mat');
DC_data_Temp_30 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_30.mat');
DC_data_Temp_32 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_32.mat');
DC_data_Temp_34 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_34.mat');
DC_data_Temp_36 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_36.mat');
DC_data_Temp_38 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_38.mat');
DC_data_Temp_40 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_40.mat');
DC_data_Temp_42 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_42.mat');
DC_data_Temp_44 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_44.mat');
DC_data_Temp_46 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_46.mat');
DC_data_Temp_48 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_48.mat');
DC_data_Temp_50 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_50.mat');
DC_data_Temp_52 = load('projects/axon_simulations/DC_temp_data/DC_data_Temp_52.mat');

% picking time shots
time1 = 5.1; % in ms
time2 = 8; % in ms
time3 = 8.5; % in ms
time4 = 9; % in ms
time5 = 9.5; % in ms
time6 = 10; % in ms
time7 = 10.5; % in ms

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
p = 0.01;

% creating a set of data from multiple experiments used to plot animation
% (the first element in the set_of_data is darkred, then the proceeding elements get
% brighter and brighter until the last element which is the brightest red)
set_of_data1 = {HH_data_Temp_default, HH_data_Temp_7, HH_data_Temp_8, HH_data_Temp_9, HH_data_Temp_10, HH_data_Temp_15, HH_data_Temp_20, HH_data_Temp_25, HH_data_Temp_30, HH_data_Temp_35};
set_of_data2 = {HH_data_Temp_30, HH_data_Temp_31, HH_data_Temp_32, HH_data_Temp_33, HH_data_Temp_34, HH_data_Temp_35};
set_of_data3 = {SC_data_Sv_200, DC_data_Sv_200};
set_of_data4 = {SC_data_Temp_20, SC_data_Temp_22, SC_data_Temp_24, SC_data_Temp_26, SC_data_Temp_28, SC_data_Temp_30, SC_data_Temp_32, SC_data_Temp_34, SC_data_Temp_36, SC_data_Temp_38, SC_data_Temp_40, SC_data_Temp_42, SC_data_Temp_44, SC_data_Temp_46, SC_data_Temp_48, SC_data_Temp_50, SC_data_Temp_52}
set_of_data5 = {DC_data_Temp_20, DC_data_Temp_22, DC_data_Temp_24, DC_data_Temp_26, DC_data_Temp_28, DC_data_Temp_30, DC_data_Temp_32, DC_data_Temp_34, DC_data_Temp_36, DC_data_Temp_38, DC_data_Temp_40, DC_data_Temp_42, DC_data_Temp_44, DC_data_Temp_46, DC_data_Temp_48, DC_data_Temp_50, DC_data_Temp_52}



% plot_animation_voltage_vs_time(DC_data, p);
% plot_animation_voltage_vs_space(DC_data, p);
% plot_animation_probabilities_vs_time(HH_data_Temp_33, p);
% plot_animation_probabilities_vs_space(HH_data_Temp_32, p);
% plot_time_and_space_shots(HH_data2, list_of_positions, list_of_times);
% plot_voltage_vs_time_comparison(set_of_data3, p);
plot_voltage_vs_space_comparison(set_of_data3, p);
% plot_Vm_and_Vm_minus_Vmy_vs_space(SC_data1, p)



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
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data, 'Vmy_all') && isfield(data, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    % plotting animation of Vm (or if SC or DC, Vm, Vmy and Vm-Vmy)
    for j=1:length(voltages)
        figure(1);
        hold on;

        xmin = 0;
        xmax = L;
        ymin = -90;
        ymax = 90;
        
        if strcmp(voltages{j}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{j}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end
            
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Axon length in cm');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:n
            x = data.(voltages{j})(i,:);
            
            % Plot the vector
            plot(t, x, 'b-');

            text(xmin + 0.05, ymax + 0.8, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
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
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data, 'Vmy_all') && isfield(data, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    for j=1:length(voltages)
        figure(1);
        hold on;

        xmin = 0;
        xmax = T;
        ymin = -90;
        ymax = 90;

        if strcmp(voltages{j}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{j}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end

        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Time in milliseconds');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:m
            x = data.(voltages{j})(:,i); 

            % Plot the vector
            plot(t, x, 'b-');
            hold on

            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

% Animation that plots the probabilities vs time
function plot_animation_probabilities_vs_space(data, p)
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 

    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = 0;
    ymax = 1;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Axon length in cm.');
    ylabel('Probabilities of gates begin open.');
    
    % Loop through each vector and plot them one by one
    for i = 1:m
        x1 = data.N_all(i,:);  
        x2 = data.M_all(i,:);  
        x3 = data.H_all(i,:);  

        % Plot the vector
        p1 = plot(t, x1, 'b-');
        hold on
        p2 = plot(t, x2, 'r-');
        hold on
        p3 = plot(t, x3, 'g-');
        hold on
        
        % Add the legend (NOTE: legends are what slows down animation)
        legend([p1, p2, p3], 'K activation gate', 'Na activation gate', 'Na inactivation date', 'Location', 'northeast');

        text(xmin, ymax, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
        
        % Add a pause to create animation effect
        pause(p);

        cla;
    end
end

% Animation that plots the probabilities vs time
function plot_animation_probabilities_vs_time(data, p)
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
    
    xmin = 0;
    xmax = T;
    ymin = -70;
    ymax = 50;
    axis([xmin xmax ymin ymax]);  % Set axis limits
    
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
    T = data_set{1}.T;
    m = data_set{1}.m;
    n = data_set{1}.n;
    dx = data_set{1}.dx;
    
    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n);
    
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data_set{1}, 'Vmy_all') && isfield(data_set{1}, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 90;
    
    for k=1:length(voltages)
        
        if strcmp(voltages{k}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{k}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end

        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Time in milliseconds');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')
        
        % Loop through each vector and plot them one by one
        for i = 1:m

            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(:,i), 'Color', [1, 0, 0] * j/length(data_set));
                hold on
            end
            
            legend('T=20C', 'T=22C', 'T=24C', 'T=26C', 'T=28C', 'T=30C', 'T=32C', 'T=34C', 'T=36C', 'T=38C', 'T=40C', 'T=42C', 'T=44C', 'T=46C', 'T=48C', 'T=50C', 'T=52C', 'Location', 'northeast');
            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

function plot_voltage_vs_space_comparison(data_set, p)
    
    % defining spatial and time variables from one of the data sets
    % (doesn't matter which since we assume they must be the same for all)
    L = data_set{1}.L;
    m = data_set{1}.m;
    n = data_set{1}.n;
    dt = data_set{1}.dt;

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, L, m); 
       
    % creating a list of voltage types to plot (either just Vm for HH model
    % or Vm, Vmy and Vm-Vmy for SC or DC model)
    voltages = {'Vm_all'};
    
    if isfield(data_set{1}, 'Vmy_all') && isfield(data_set{1}, 'Vm_minus_Vmy')
        voltages{end+1} = 'Vmy_all';
        voltages{end+1} = 'Vm_minus_Vmy';
    end
    
    figure(1);
    hold on;

    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    for k=1:length(voltages)
        
        if strcmp(voltages{k}, 'Vm_all')
            voltage_name = 'V_m';
        elseif strcmp(voltages{k}, 'Vmy_all')
            voltage_name = 'V_{my}';
        else
            voltage_name = 'V_m - V_{my}';
        end
    
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Length of axon in cm');
        ylabel(['$', voltage_name, '$ in millivolts.'], 'Interpreter', 'latex')

        % Loop through each vector and plot them one by one
        for i = 1:n
            for j = 1:length(data_set)
                plot(t, data_set{j}.(voltages{k})(i,:), 'Color', [1, 0, 0] * j/length(data_set));
                hold on
            end

            % Add the legend (NOTE: the legend is what is slowing down the animation)
            legend('T=20C', 'T=22C', 'T=24C', 'T=26C', 'T=28C', 'T=30C', 'T=32C', 'T=34C', 'T=36C', 'T=38C', 'T=40C', 'T=42C', 'T=44C', 'T=46C', 'T=48C', 'T=50C', 'T=52C', 'Location', 'northeast');
            text(xmin, ymax + 0.1, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);

            cla;
        end
    end
end

% plotting Vm vs Space and Vm-Vmy vs Space in the same plot (only for SC and DC)
function plot_Vm_and_Vm_minus_Vmy_vs_space(data, p)
    
    % defining spatial and time variables from one of the data sets
    L = data.L;
    m = data.m;
    n = data.n;
    dt = data.dt;

    % x axis is the axon time
    t = linspace(0, L, m); 
    
    figure(1);
    
    xmin = 0;
    xmax = L;
    ymin = -90;
    ymax = 90;
    
    % Loop through each vector and plot them one by one
    for i = 1:n
        cla;
        
        plot(t, data.Vm_all(i,:), 'b-');
        hold on
        plot(t, data.Vm_minus_Vmy(i,:), 'r-');
        hold off
        
        axis([xmin xmax ymin ymax]);  % Set axis limits
        xlabel('Length of axon in cm');
        ylabel('Voltage in millivolts.', 'Interpreter', 'latex')
        
        % Add the legend (NOTE: the legend is what is slowing down the animation)
        % legend([p1, p2], 'data1 voltage', 'data2 voltage', 'Location', 'northeast');
        text(xmin+0.05*(xmax-xmin), ymax-5, sprintf('Time: %.3f ms', round(i*dt, 3)), 'FontSize', 12, 'BackgroundColor', 'w');
    
        % Add a pause to create animation effect
        pause(p);

    end
end


% Plotter function used to plot space shots of multiple axon segments of 
% data vs space where data could be one of the following: Vm, Vmy, 
% Vm_minus_Vmy, n, m, h
% Kevin Roberts
% May 2025

function plot_data_vs_space_axon_segments(data, data_type, axon_segments, time_shot)
    
    % Unpack values
    L = data.L;         % Total axon length (cm)
    m = data.m;         % Number of spatial points
    dt = data.dt;       % Time step (ms)
    dx = data.dx;       % Spatial step (cm)
    
    x_axis = 'Length of axon in cm';

    [data_types, display_names, y_axis, ymin, ymax] = helper_functions.get_plot_names_and_limits(data_type);
    data_type = data_types{1}; % there's only one for this plotter function
    display_name = display_names{1}; % again, only one here

    figure(1)
    hold on


    interval1 = axon_segments{1};
    interval2 = axon_segments{2};
    interval3 = axon_segments{3};

    % Time vector
    total_timesteps = size(data.(data_type), 1);
    t_vec = (0:total_timesteps-1) * dt;

    % --- Interval 1: Axon Segment 5 ---
    x_positions1 = linspace(interval1(1), interval1(2), 4);
    spatial_indices1 = round(x_positions1 / dx) + 1;

    % Interval range for legend
    interval1_range = sprintf(['%.2f-%.4f cm of ', display_name], interval1(1), interval1(2));

    for i = 1:length(spatial_indices1)
        xi = spatial_indices1(i);
        Vm_trace = data.(data_type)(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'b-', 'DisplayName', ['Axon Segment 5: ', interval1_range]);
        else
            plot(t_vec, Vm_trace, 'b-', 'HandleVisibility', 'off');
        end
    end

    % --- Interval 2: Axon Segment 10 ---
    x_positions2 = linspace(interval2(1), interval2(2), 4);
    spatial_indices2 = round(x_positions2 / dx) + 1;

    % Interval range for legend
    interval2_range = sprintf(['%.2f-%.4f cm of ', display_name], interval2(1), interval2(2));

    for i = 1:length(spatial_indices2)
        xi = spatial_indices2(i);
        Vm_trace = data.(data_type)(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'r-', 'DisplayName', ['Axon Segment 10: ', interval2_range]);
        else
            plot(t_vec, Vm_trace, 'r-', 'HandleVisibility', 'off');
        end
    end

    % --- Interval 3: Axon Segment 15 ---
    x_positions3 = linspace(interval3(1), interval3(2), 4);
    spatial_indices3 = round(x_positions3 / dx) + 1;

    % Interval range for legend
    interval3_range = sprintf(['%.2f-%.4f cm of ', display_name], interval3(1), interval3(2));

    for i = 1:length(spatial_indices3)
        xi = spatial_indices3(i);
        Vm_trace = data.(data_type)(:, xi);
        if i == 1
            plot(t_vec, Vm_trace, 'k-', 'DisplayName', ['Axon Segment 15: ', interval3_range]);
        else
            plot(t_vec, Vm_trace, 'k-', 'HandleVisibility', 'off');
        end
    end

    legend('show', 'Location', 'northeast')

    hold off
    
    xlabel(x_axis, 'Interpreter', 'latex');
    ylabel(y_axis, 'Interpreter', 'latex');
    
end 
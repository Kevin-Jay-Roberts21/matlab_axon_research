% Plotter function used to plot one time shot of data vs space where 
% data could be one of the following: Vm, Vmy, Vm_minus_Vmy, n, m, h, nmh,
% Vm_and_Vm_minus_Vmy, or Vm_and_Vmy_and_Vm_minus_Vmy
% Kevin Roberts
% May 2025

function plot_data_vs_space_at_one_time_shot(data, data_type, time_shot)

    L = data.L;
    m = data.m;
    dt = data.dt;

    t = linspace(0, L, m); 
    x_axis = 'Length of axon in cm';
    xmin = 0;
    xmax = L;

    [data_types, display_names, y_axis, ymin, ymax] = helper_functions.get_plot_names_and_limits(data_type);

    figure(1);
    clf;
    hold on;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel(x_axis, 'Interpreter', 'latex');
    ylabel(y_axis, 'Interpreter', 'latex');

    % Initial plots to get handles for each line
    plot_handles = gobjects(1, length(data_types));
    for j = 1:length(data_types)
        x = data.(data_types{j})(round(time_shot/dt),:); 
        plot_handles(j) = plot(t, x, 'DisplayName', sprintf([display_names{j}, ' at t = %.2f ms'], time_shot));
    end

    % Show the legend once
    legend('Interpreter', 'latex');

end 
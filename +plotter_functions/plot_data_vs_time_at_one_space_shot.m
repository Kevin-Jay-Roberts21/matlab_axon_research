% Plotter function used to plot one space shot of data vs time where 
% data could be one of the following: Vm, Vmy, Vm_minus_Vmy, n, m, h, nmh,
% Vm_and_Vm_minus_Vmy, or Vm_and_Vmy_and_Vm_minus_Vmy
% Kevin Roberts
% May 2025

function plot_data_vs_time_at_one_space_shot(data, data_type, space_shot)

    T = data.T;
    n = data.n;
    dx = data.dx;

    t = linspace(0, T, n); 
    x_axis = 'Time in milliseconds';
    xmin = 0;
    xmax = T;

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
        x = data.(data_types{j})(:,round(space_shot/dx)); 
        plot_handles(j) = plot(t, x, 'DisplayName', sprintf([display_names{j}, ' at x = %.4f cm'], space_shot));
    end

    % Show the legend once
    legend('Interpreter', 'latex');
   

end 
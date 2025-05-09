% Plotter function used to plot the animation of the data vs time where
% data could be one of the following: Vm, Vmy, Vm_minus_Vmy, n, m, h,
% n_and_m_and_h, Vm_and_Vm_minus_Vmy, or Vm_and_Vmy_and_Vm_minus_Vmy
% Kevin Roberts
% May 2025

function plot_animation_data_vs_time(data, data_type, p)

    T = data.T;
    m = data.m;
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
        x = data.(data_types{j})(:,1); 
        plot_handles(j) = plot(t, x, 'DisplayName', display_names{j});
    end

    % Show the legend once
    legend('Interpreter', 'latex');

    % Create and store the text handle
    text_handle = text(xmin + 0.01*(xmax - xmin), ymax - 0.05*(ymax - ymin), ...
        '', 'FontSize', 12, 'BackgroundColor', 'w', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % Animation loop — update data instead of re-plotting
    for i = 1:m
        for j = 1:length(data_types)
            plot_handles(j).YData = data.(data_types{j})(:,i);
        end

        text_handle.String = sprintf('Space: %.5f cm', round(i*dx, 5));
        pause(p);
    end
end
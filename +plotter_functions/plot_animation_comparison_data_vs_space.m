% Plotter function used to plot the animation comparison of multiple lines 
% of data vs space where data could be one of the following: Vm, Vmy, 
% Vm_minus_Vmy, n, m, h.
% Kevin Roberts
% May 2025

function plot_animation_comparison_data_vs_space(data_set, data_set_names, data_type, p)

    L = data_set{1}.L;
    m = data_set{1}.m;
    dt = data_set{1}.dt;

    t = linspace(0, L, m); 
    x_axis = 'Length of axon in cm';
    xmin = 0;
    xmax = L;

    [data_types, display_names, y_axis, ymin, ymax] = helper_functions.get_plot_names_and_limits(data_type);
    data_type = data_types{1}; % there's only one for this plotter function
    display_name = display_names{1}; % again, only one here
    
    figure(1);
    clf;
    hold on;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel(x_axis, 'Interpreter', 'latex');
    ylabel(y_axis, 'Interpreter', 'latex');

    % Initial plots to get handles for each line
    plot_handles = gobjects(1, length(data_set));
    for j = 1:length(data_set)
        x = data_set{j}.(data_type)(1,:); 
        safe_name = helper_functions.escape_latex_chars(data_set_names{j});
        full_name = [display_name, ' of ', safe_name];
        plot_handles(j) = plot(t, x, 'DisplayName', full_name);
    end

    % Show the legend once
    legend(plot_handles, 'Interpreter', 'latex');

    % Create and store the text handle
    text_handle = text(xmin + 0.01*(xmax - xmin), ymax - 0.05*(ymax - ymin), ...
        '', 'FontSize', 12, 'BackgroundColor', 'w', 'Interpreter', 'latex', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % Animation loop — update data instead of re-plotting
    for i = 1:m
        for j = 1:length(data_set)
            plot_handles(j).YData = data_set{j}.(data_type)(i,:);
        end

        text_handle.String = sprintf('Space: %.2f cm', round(i*dt, 2));
        pause(p);
    end

end 
% Plotter function used to plot multiple space shots of data vs time where 
% data could be one of the following: Vm, Vmy, Vm_minus_Vmy, n, m, h.
% Kevin Roberts
% May 2025

function plot_data_vs_time_at_space_shots(data, data_type, space_shots)

    T = data.T;
    n = data.n;
    dx = data.dx;
    
    t = linspace(0, T, n); 
    x_axis = 'Time in milliseconds';
    xmin = 0;
    xmax = T;

    [data_types, display_names, y_axis, ymin, ymax] = helper_functions.get_plot_names_and_limits(data_type);
    data_type = data_types{1}; % there's only one for this plotter function
    display_name = display_names{1}; % again, only one here

    figure(1);
    clf;
    hold on;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel(x_axis, 'Interpreter', 'latex');
    ylabel(y_axis, 'Interpreter', 'latex');

    for i = 1:length(space_shots)
        plot(t, data.(data_type)(:,round(space_shots{i}/dx)))
        hold on
    end
    
    % describing plots using legends
    legendStrings1 = {};
    for i  = 1:length(space_shots)
        legendStrings1{end+1} = sprintf([display_name, ' at x = %.4f cm'], space_shots{i});
    end
    legend(legendStrings1, 'Interpreter','latex', 'FontSize', 12)

end 
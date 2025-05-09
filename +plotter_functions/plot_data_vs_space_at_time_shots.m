% Plotter function used to plot multiple time shots of data vs space where 
% data could be one of the following: Vm, Vmy, Vm_minus_Vmy, n, m, h.
% Kevin Roberts
% May 2025

function plot_data_vs_space_at_time_shots(data, data_type, time_shots)

    L = data.L;
    m = data.m;
    dt = data.dt;
    
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

    for i = 1:length(time_shots)
        plot(t, data.(data_type)(round(time_shots{i}/dt),:))
        hold on
    end
    
    % describing plots using legends
    legendStrings1 = {};
    for i  = 1:length(time_shots)
        legendStrings1{end+1} = sprintf([display_name, ' at t = %.2f ms'], time_shots{i});
    end
    legend(legendStrings1, 'Interpreter','latex', 'FontSize', 12)

end 
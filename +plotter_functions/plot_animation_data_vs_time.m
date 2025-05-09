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

    % TEMPORAL PROFILE %
    % x axis is the axon time
    t = linspace(0, T, n); 
    
    xmin = 0;
    xmax = T;
    ymin = -90;
    ymax = 90;
    
    data_types = {};
    display_names = {};
    y_axis = '';
    
    if strcmp(data_type, "Vm")
        data_types = {'Vm_all'};
        display_names = {'$V_m$'};
        y_axis = 'Voltage in mV';
    elseif strcmp(data_type, "Vmy")
        data_types = {'Vmy_all'};
        display_names = {'$V_{my}$'};
        y_axis = 'Voltage in mV';
    elseif strcmp(data_type, "Vm_minus_Vmy")
        data_types = {'Vm_minus_Vmy'};
        display_names = {'$V_m - V_{my}$'};
        y_axis = 'Voltage in mV';
    elseif strcmp(data_type, "n")
        data_types = {'N_all'};
        display_names = {'$n$'};
        y_axis = 'Gating variable $n$';
        ymin = 0;
        ymax = 1;
    elseif strcmp(data_type, "m")
        data_types = {'M_all'};
        display_names = {'$m$'};
        y_axis = 'Gating variable $m$';
        ymin = 0;
        ymax = 1;
    elseif strcmp(data_type, "H")
        data_types = {'H_all'};
        display_names = {'$h$'};
        y_axis = 'Gating variable $h$';
        ymin = 0;
        ymax = 1;
    elseif strcmp(data_type, "nmh")
        data_types = {'N_all', 'M_all', 'H_all'};
        display_names = {'$n$', '$m$', '$h$'};
        y_axis = 'Gating variables $n$, $m$, and $h$';
        ymin = 0;
        ymax = 1;
    elseif strcmp(data_type, "Vm_and_Vm_minus_Vmy")
        data_types = {'Vm_all', 'Vm_minus_Vmy'};
        display_names = {'$V_m$', '$V_m - V_{my}$'};
        y_axis = 'Voltage in mV';
    elseif strcmp(data_type, "Vm_and_Vmy_and_Vm_minus_Vmy")
        data_types = {'Vm_all', 'Vmy_all', 'Vm_minus_Vmy'};
        display_names = {'$V_m$', '$V_{my}$', '$V_m - V_{my}$'};
        y_axis = 'Voltage in mV';
    else
        disp('You did not pick a valid data_type name. Try again.');
        return  % exits the script
    end
    
    
    figure(1);
    hold on;

    axis([xmin xmax ymin ymax]);  % Set axis limits
    xlabel('Time in milliseconds');
    ylabel(y_axis, 'Interpreter', 'latex')

    for i = 1:m
        
        for j=1:length(data_types)
            x = data.(data_types{j})(:,i); 
            
            plot(t, x, 'b-');
            hold on

            text(xmin + 0.2, ymax + 0.1, sprintf('Space: %.5f cm', round(i*dx, 5)), 'FontSize', 12, 'BackgroundColor', 'w');

            % Add a pause to create animation effect
            pause(p);
            
        end
        
        cla;
    end

end 
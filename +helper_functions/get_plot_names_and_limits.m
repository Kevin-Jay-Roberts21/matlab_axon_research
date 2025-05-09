% Function created to return names of axis, data, axis limits, etc.
% Kevin Roberts
% May 2025

function [data_types, display_names, y_axis, ymin, ymax] = get_plot_names_and_limits(data_type)
    
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
    
end
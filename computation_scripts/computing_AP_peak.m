% Code used to compute the peak voltage of an action potential for HH and
% SC/DC models
% Kevin Roberts
% July 2025
function list_of_ap_peaks = computing_AP_peak(data)
    
    list_of_ap_peaks = zeros(1, length(data));

    for i = 1:length(data)
    
        % SCHEME
        Vm_all = data{i}.Vm_all;
        m = data{i}.m;
        
        % Choose the spatial midpoint
        mid_index = round(m/2);
        
        % Extract voltage at the midpoint over time
        Vm_time_at_mid = Vm_all(:, mid_index);
        
        % Equilibrium voltage is assumed to be the initial voltage
        V_eq = Vm_time_at_mid(1);
        
        % Find the peak voltage and compute half-height voltage
        [V_peak, peak_index] = max(Vm_time_at_mid);

        list_of_ap_peaks(i) = V_peak;
    end
end
% Code used to Computed the duration of an action potential for HH and
% SC/DC models
% Kevin Roberts
% July 2025

function list_of_ap_durations = computing_AP_duration(data)
    
    list_of_ap_durations = zeros(1, length(data));

    for i = 1:length(data)
        Vm_all = data{i}.Vm_all;
        dt = data{i}.dt;
        m = data{i}.m;
        
        % Choose the spatial midpoint
        mid_index = round(m/2);
        
        % Extract voltage at the midpoint over time
        Vm_time_at_mid = Vm_all(:, mid_index);
        
        % Equilibrium voltage is assumed to be the initial voltage
        V_eq = Vm_time_at_mid(1);
        
        % Find the peak voltage and compute half-height voltage
        [V_peak, peak_index] = max(Vm_time_at_mid);
        V_half = V_eq + 0.5 * (V_peak - V_eq);
        
        % Find the first index where voltage is closest to V_half before the peak
        [~, t_h1_index] = min(abs(Vm_time_at_mid(1:peak_index) - V_half));
        
        % Find the second index after the peak where voltage is again closest to V_half
        [~, idx2_rel] = min(abs(Vm_time_at_mid(peak_index:end) - V_half));
        t_h2_index = peak_index - 1 + idx2_rel;  % convert relative index to absolute index
        
        % Convert to actual time values
        ap_duration = 2*dt*(t_h2_index - t_h1_index);

        list_of_ap_durations(i) = ap_duration;
    end

end

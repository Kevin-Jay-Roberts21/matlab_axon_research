function [speed,time_difference, voltage_difference] = repolarization(Uall, position, equilibrium, h, k)
    
    
    vector_of_voltages = Uall(:,round(position/h));
    max_voltage = max(vector_of_voltages);
    [~, time_of_max_voltage] = min(abs(vector_of_voltages - max_voltage));
    final_time = 0;
    for i = time_of_max_voltage:length(vector_of_voltages)
        voltage_value = vector_of_voltages(i);
        if voltage_value < equilibrium
            final_time = (i-1)*k;
            break
        end
    end
    time_of_max_voltage = time_of_max_voltage*k % converting to ms
    final_time
    equilibrium
    max_voltage
    time_difference = abs(time_of_max_voltage - final_time);
    voltage_difference = abs(max_voltage - equilibrium);
    speed = voltage_difference/time_difference;
end
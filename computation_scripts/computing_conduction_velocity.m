% A MATLAB script to compute the conduction velocity of the voltage
% simulations. More accurate results for longer simulations.
% Kevin Roberts
% February 2025

function final_cvs = computing_conduction_velocity(spaces, data)
    
    sum_of_cvs = zeros(1, length(data));

    for j = 1:size(spaces, 1)
        
        list_of_cv = zeros(1, length(data)); % Preallocate for efficiency
        
        for i = 1:length(data)
            
            x1 = spaces(j, 1);
            x2 = spaces(j, 2);
            
            % need to find the index of x1 and x2
            index_x1 = x1/data{i}.dx;
            index_x2 = x2/data{i}.dx;
            
            % identifying the space vectors at index_t1 and index_t2
            vec_t1 = data{i}.Vm_all(:,index_x1);
            vec_t2 = data{i}.Vm_all(:,index_x2);
        
            % computing the positions of where the max voltage is at t1 and t2:
            [Vm_at_x1, index_t1] = max(vec_t1); % Time index where voltage peaks at x1
            [Vm_at_x2, index_t2] = max(vec_t2); % Time index where voltage peaks at x2
            
            % Note that x1 and x2 are just the indices. We need to find their actual
            % spatial position in cm. This is done by identifiying the mesh from the
            % data
            t1 = index_t1*data{i}.dt; % (in cm)
            t2 = index_t2*data{i}.dt; % (in cm)
            
            % finally, compute the conduction velocity
            cv = (x2 - x1)/(t2 - t1) * 10; % *10 to convert to m/s
    
            list_of_cv(i) = cv;
        end
        
        sum_of_cvs = sum_of_cvs + list_of_cv;
    end

    % now dividing sum_of_cvs by the number of cv we computed
    final_cvs = sum_of_cvs/size(spaces, 1);

end




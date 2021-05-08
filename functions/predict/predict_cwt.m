function [predictions] = predict_cwt(velocity, wavelets, params)
%predict_cwt
%
% Syntax: predictions = predict_cwt(velocity, wavelets, params)
%
% Implements cwt predictions.
narginchk(3, 3)
nargoutchk(1, 1)
velocity.time = 0;
validatevelocity(velocity)

% Convienence variables
n_samples = size(velocity.xt, 1);
n_wavelets = size(wavelets.even, 1);
n_sensors = params.sensors.n_sensors;
predictions = zeros(n_samples, 3);

% Apply input mode
[vx, vy] = cwt_input_mode_velocity(velocity, params);
vx = permute(vx, [1 3 2]);
vy = permute(vy, [1 3 2]);
clear velocity;

% Compute coefficients
even = wavelets.even; 
odd = wavelets.odd; 
navelet = wavelets.navelet;
X = wavelets.sources(:, 1);
Y = wavelets.sources(:, 2);
scaling = 1 ./ sqrt(abs(params.sensors.y_value - Y))';
clear wavelets;

% Compute sizes (8 sensors)
% How much memory does each worker use?
worker_copy_size = ...
    4 * n_wavelets +... % magnitude, X, Y, filter
    5 * n_sensors +...  % vx, vy, wavelets
    8;                  % prediction, coefs
% How much memory can the CWT use?
n_workers = params.memory.max_workers;
system_buffer = 2; % Reserve 2G for the system
matlab_buffer = n_workers * 0.95 + 2; % Reserve 0.95G for each worker and 2G for the interface
available_memory = params.memory.mem - system_buffer - matlab_buffer;

% How much memory remains for the iteration?
iteration_memory = available_memory - (n_workers + 1) * (worker_copy_size * 8 / 1024^3);

one_iteration_size = ...
    4 * n_wavelets + ... % x_odd, x_even, y_nav, y_odd
    2 * n_sensors;       % vx, vy

% Determine how many sources can be computed in one iteration
factor = 1.8;
step_size = floor(iteration_memory / (factor * (one_iteration_size * 8 / 1024^3)));
starts = 1:step_size:n_samples;          
%               
% estimated_memory_usage = ...
%     system_buffer + ...
%     matlab_buffer + (n_workers + 1) * (worker_copy_size * 8 / 1024^3) + ...
%     min(step_size, n_samples) * factor * (one_iteration_size * 8 / 1024^3);
% disp(['Using ' num2str(step_size) ' locations per iteration'])
% disp(['    Iterations: ' num2str(length(starts))]) 
% disp(['    Estimated memory: ' num2str(estimated_memory_usage) ' GB'])

% Compute predictions in iterations
for idx = starts
    indices = idx:min(idx+(step_size-1), n_samples);
    n_indices = length(indices);
    
    % Compute coefficients (at training points)
    x_even = scaling .* (vx(indices,:) * even');
    x_odd =  scaling .* (vx(indices,:) * odd'); 
    y_nav =  scaling .* (vy(indices,:) * navelet');
    y_odd =  scaling .* (vy(indices,:) * odd');
    
    % Combine all coefficients
    magnitude = sqrt(x_even.^2 + y_nav.^2 + x_odd.^2 + y_odd.^2);
    magnitude = magnitude ./ max(magnitude, [], 2);
    
    loop_predictions = zeros(n_indices, 3);
    parfor index_idx = 1:n_indices
        % Allow matlab to slice the matrix, limits memory usage
        par_mag = magnitude(index_idx,:);
        par_vx = vx(index_idx, :);
        par_vy = vy(index_idx, :);
        
        filter = par_mag > params.cwt.threshold_min & par_mag < params.cwt.threshold_max;
        n_sources = sum(filter(:));
        if n_sources <= 0
            % warning('CWT:filterExhausted', ...
            %     ['No coefficient values passed the threshold filter. '...
            %      'Used the maximum value instead.'])
            [~, i] = max(par_mag(:));
            p_hat = [X(i), Y(i)];
        elseif n_sources < 7 % The number of variables that is fitted
%             warning('CWT:filterTooRestrictive', ...
%                 ['Not enough points passed the threshold filter to fit a Gaussian.\n',...
%                  'Used the maximum value within the filter instead.'])
            [~, i] = max(par_mag(filter));
            X_filter = X(filter);
            Y_filter = Y(filter);
            p_hat = [X_filter(i), Y_filter(i)];
        else
            res = fmgaussfit(X(filter), Y(filter), par_mag(filter)', params);
            p_hat = res(5:6);
        end
        
        % Compute orientation
        p_hat_wavelets = compute_wavelets(p_hat, params);
        p_hat_scaling = 1 ./ sqrt(abs(params.sensors.y_value - p_hat_wavelets.sources(:,2)));
        p_hat_x_even = p_hat_scaling .* sum(par_vx .* p_hat_wavelets.even, 2);
        p_hat_x_odd =  p_hat_scaling .* sum(par_vx .* p_hat_wavelets.odd, 2);
        p_hat_y_nav =  p_hat_scaling .* sum(par_vy .* p_hat_wavelets.navelet, 2);
        p_hat_y_odd =  p_hat_scaling .* sum(par_vy .* p_hat_wavelets.odd, 2);
        
        phi_x = mod(atan2(params.cwt.c_x * p_hat_x_odd, p_hat_x_even), 2*pi); 
        phi_y = mod(atan2(params.cwt.c_y * p_hat_y_nav, p_hat_y_odd), 2*pi);

        switch params.sensors.input_mode
            case 'x'
                loop_predictions(index_idx, :) = [p_hat phi_x];
            case 'y'
                loop_predictions(index_idx, :) = [p_hat phi_y];
            otherwise % x+y and x|y  
                loop_predictions(index_idx, :) = [p_hat compute_angle_average([phi_x, phi_y])];
        end
    end
    predictions(indices,:) = loop_predictions;
end

end
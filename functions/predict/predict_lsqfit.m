function [predictions] = predict_lsqfit(velocity, params)
%predict_lsqfit
%
% Syntax: predictions = predict_lsqfit(velocity, params)
%
% Implements lsqfit predictions.
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(velocity)

% Easy access in parfor loop (reduce memcopy)
n_samples = size(velocity.input, 1);
x_range = params.domain.x_range;
y_range = params.domain.y_range;
a_range = params.domain.azimuth_range;
input_velocity = velocity.input;
sensor_locations = params.sensors.locations.x;
lsq = params.lsq;
% Preallocate results
predictions  = zeros(n_samples, 3);

% Do parallel work
parfor sample_idx = 1:n_samples
    % Extract iteration data (reduce memcopy)
    sample_velocity = input_velocity(sample_idx, :); % is a row vector!
       
    x0 = lsq.starting_estimate;
    
    fit_function = @(x, xdata) simulate_lsq_velocity(x, params);
    options = optimoptions('lsqcurvefit',...
        'MaxIterations', lsq.max_iterations,...
        'FiniteDifferenceStepSize', lsq.gradient_step,...
        'FiniteDifferenceType', 'central',...
        'StepTolerance', lsq.step_tolerance,...
        'FunctionTolerance', lsq.function_tolerance,...
        'OptimalityTolerance', lsq.optimality_tolerance,...
        'Display', 'off');
    x_fit = lsqcurvefit(fit_function, x0, sensor_locations, sample_velocity,...
        [min(x_range), min(y_range), min(a_range)],... lower bound
        [max(x_range), max(y_range), max(a_range)],... upper bound
        options);
    
    % Write predictions
    predictions(sample_idx, :) = x_fit;
end

end
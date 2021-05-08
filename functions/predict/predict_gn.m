function [predictions, estimates] = predict_gn(velocity, params)
%predict_gn
%
% Syntax: predictions = predict_gn(velocity, params)
%
% Implements GN predictions.
narginchk(2, 2)
nargoutchk(1, 2)
validatevelocity(velocity)

fun = @(input, estimates) predict(input, estimates, params);
[predictions, estimates] = compute_nonlin_prediction(velocity, fun, params.gn.max_iterations);

predictions = determine_orientation_phase(predictions, velocity, params);

end

function [prediction, estimates] = predict(input, estimates, params)
%predict
%
% Syntax:  [prediction, estimates] = gn_predict(input, estimates, params)
%
% Performs a single prediction using the GN algorithm.
%
% params.gn contains the method's hyperparameters:
% - params.gn.starting_estimate    is the starting estimate [x y azimuth]
% - params.gn.gradient_step        is the step used in computing the derivatives
%                                  numerically
% - params.gn.max_iterations       is a stopping condition on the number of
%                                  iterations
% - params.gn.min_step_tolerance   is a stopping condition on the change in estimate
% - params.gn.step_size            is the minimum step size in the update rule
% - params.gn.norm_limit           is the maximum norm of an update step.

narginchk(3, 3)
nargoutchk(1, 3)

gn = params.gn;
validateattributes(input, {'double'}, {'ncols', 1}, 'gn_predict', 'input', 1)
validateattributes(estimates, {'double'}, {'2d', 'nrows', gn.max_iterations, 'ncols', 3}, 'gn_predict', 'estimates', 2)
estimates(1,:) = gn.starting_estimate;

% Run the iterations
lastwarn('');
for idx = 2:gn.max_iterations
    col = idx - 1;
    
    estimate = estimates(col, :);
    
    % Compute velocity and gradient at the estimate
    velocity = compute_sensor_velocity(estimate, 0, params);
    velocity = apply_input_mode(velocity, params);
    grad = compute_gradient_numerical(estimate, params, gn.gradient_step);
    direction = (grad' * grad) \ grad' * (input - squeeze(abs(velocity.input)));
    [~, id] = lastwarn;
    if strcmp(id, 'MATLAB:nearlySingularMatrix')
        estimates(idx,:) = estimates(col,:);
        break;
    end
    estimates(idx,:) = update_nonlin_estimate(estimate, direction, gn.step_size, gn.norm_limit, params);
    
    % Check stop criteria
    if norm(estimates(idx,:) - estimates(col,:)) < gn.step_tolerance
        break
    end
end

prediction = estimates(idx,:);

end
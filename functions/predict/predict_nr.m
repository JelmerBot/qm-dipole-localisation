function [predictions, estimates] = predict_nr(velocity, params)
%predict_gn
%
% Syntax: predictions = predict_nr(velocity, params)
%
% Implements NR predictions.
narginchk(2, 2)
nargoutchk(1, 2)
validatevelocity(velocity)

fun = @(input, estimates) predict(input, estimates, params);
[predictions, estimates] = compute_nonlin_prediction(velocity, fun, params.nr.max_iterations);

predictions = determine_orientation_phase(predictions, velocity, params);

end

function [prediction, estimates] = predict(input, estimates, params)
%predict
%
% Syntax:  [prediction, estimates] = gn_predict(input, estimates, params)
%
% Performs a single prediction using the NR algorithm.
%
% params.nr contains the method's hyperparameters:
% - params.nr.starting_estimate    is the starting estimate [x y azimuth]
% - params.nr.gradient_step        is the step used in computing the derivatives
%                                  numerically
% - params.nr.max_iterations       is a stopping condition on the number of
%                                  iterations
% - params.nr.min_step_tolerance   is a stopping condition on the change in estimate
% - params.nr.step_size            is the minimum step size in the update rule
% - params.nr.max_step             is the maximum norm of an update step.
% - params.nr.regularization       regulates the hessian (H = H + aI).
narginchk(3, 3)
nargoutchk(1, 3)

nr = params.nr;
validateattributes(input, {'double'}, {'ncols', 1}, 'gn_predict', 'input', 1)
validateattributes(estimates, {'double'}, {'2d', 'nrows', nr.max_iterations, 'ncols', 3}, 'gn_predict', 'estimates', 2)
estimates(1,:) = nr.starting_estimate;
% regu = diag(repmat(nr.regularization, 3, 1));

% negative_counter = zeros(1, nr.max_iterations);

% Run the iterations
lastwarn('');
for idx = 2:nr.max_iterations
    col = idx - 1;
    
    estimate = estimates(col, :);
    
    % Estimate g and G numerically
    g = g_numerical(input, estimate, params);
    G = G_numerical(input, estimate, params);
    
%     negative_counter(idx) = any(eig(G) < 0);
    
%     direction = (G + regu) \ g;
    direction = G \ g;
    [~, id] = lastwarn;
    if strcmp(id, 'MATLAB:nearlySingularMatrix')
        estimates(idx,:) = estimates(col,:);
        break;
    end
    % NR updates in -direction instead of +direction!!!
    estimates(idx,:) = update_nonlin_estimate(estimate, -direction, nr.step_size, nr.norm_limit, params);
    
    % Check stop criteria
    if norm(estimates(idx,:) - estimates(col,:)) < nr.step_tolerance
        break
    end
end

prediction = estimates(idx,:);
% all(negative_counter(2:idx))

end

function [g] = g_numerical(input, estimate, params)
    nr = params.nr;
    velocity = compute_sensor_velocity(estimate, 0, params);
    velocity = apply_input_mode(velocity, params);
    grad = compute_gradient_numerical(estimate, params, nr.gradient_step);
    
    g = grad' * (input - squeeze(abs(velocity.input)));
end

function [G] = G_numerical(input, estimate, params)
    nr = params.nr;
    x_values = [estimate(1) - nr.gradient_step, estimate(1), estimate(1) + nr.gradient_step];
    y_values = [estimate(2) - nr.gradient_step, estimate(2), estimate(2) + nr.gradient_step];
    azi_values = [estimate(3) - nr.gradient_step, estimate(3), estimate(3) + nr.gradient_step];
    [Y, X, A1] = ndgrid(y_values, x_values, azi_values);
    G_sources = [X(:) Y(:) A1(:)];
    n_sources = size(G_sources, 1);
    clear X Y Z A1

    g = zeros(n_sources, 3);
        
    for row = 1:n_sources
        g(row, :) = g_numerical(input, G_sources(row, :), params);
    end
    
    x = reshape(permute(g, [1 3 2]), [3 3 3  3]);
    [dx, dy, da] = gradient(x, nr.gradient_step, nr.gradient_step, nr.gradient_step, 1);
    G = cat(2, squeeze(dx(2,2,2,:)), squeeze(dy(2,2,2,:)), squeeze(da(2,2,2,:)));
end

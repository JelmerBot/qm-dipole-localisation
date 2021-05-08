function [predictions, estimates] = predict_qm(velocity, params)
%predict_qm
%
% Syntax: [predictions, estimates] = predict_qm(velocity, params)
%
% Implements Quadrature predictions.
narginchk(2, 2)
nargoutchk(1, 2)
validatevelocity(velocity)

if ~strcmp(params.sensors.input_mode, 'x+y')
   error('Quadrature only works with vx and vy!') 
end

% Preallocation and easy access
n_samples = size(velocity.xt, 1);
velocity_xt = velocity.xt;
velocity_yt = velocity.yt;
predictions  = zeros(n_samples, 3);
estimates = zeros(n_samples, 3);

% Compute quadrature
quadrature = sqrt(velocity.xt.^2 + 0.5 .* velocity.yt.^2);
quadrature = permute(quadrature, [1 3 2]);
quadrature_norm = quadrature ./ max(quadrature, [], 2);

% Feature values
rho_anch = 2 / sqrt(5);
features.rho_value = sqrt((1 + 0.5 * rho_anch^2 + 4 * rho_anch^4) /...
                            (1 + rho_anch^2)^5);
features.distance_scaling = 2 * rho_anch;

parfor sample_idx = 1:n_samples
    sample_quadrature_norm = quadrature_norm(sample_idx, :);
    sample_quadrature = quadrature(sample_idx, :);
    sample_vx = velocity_xt(sample_idx, :);
    sample_vy = velocity_yt(sample_idx, :);
    
    estimate = qm_feature_estimate(sample_quadrature_norm, features, params.sensors.locations.x);
    estimate = qm_estimate_orientation(sample_vx, sample_vy, estimate, params);
    fit_ori = estimate;

    for repeat_idx = 1:params.qm.refine_repeats
        fit_loc = qm_fit_location(sample_quadrature, fit_ori, params);
        fit_ori = qm_estimate_orientation(sample_vx, sample_vy, fit_loc, params);
    end
    
    estimates(sample_idx,:) = estimate;
    predictions(sample_idx,:) = fit_ori;
end

end

function prediction = qm_feature_estimate(quadrature_norm, features, sensor_locations)

n_sensors = length(sensor_locations);
[~, x_peak] = max(quadrature_norm, [], 2);
quadrature_thresholded = quadrature_norm > features.rho_value;

candidates = find(~quadrature_thresholded);
left_candidates = candidates(candidates < x_peak);
if isempty(left_candidates)
    left = sensor_locations(1);
    left_intersected = false;
else
    left = left_candidates(end);
    % Finds the x-location where quadrature_norm crosses features.rho_value
    % using linear interpolation.
    left = interp1(quadrature_norm([left, left+1]),...
                    sensor_locations([left, left+1]),...
                    features.rho_value, 'linear');
    left_intersected = true;
end
right_candidates = candidates(candidates > x_peak);
if isempty(right_candidates)
    right = sensor_locations(n_sensors);
    right_intersected = false;
else
    right = right_candidates(1);
    % Finds the x-location where quadrature_norm crosses features.rho_value
    % using linear interpolation.
    right = interp1(quadrature_norm([right-1, right]),...
                    sensor_locations([right-1, right]),...
                    features.rho_value, 'linear');
    right_intersected = true;
end

b = sensor_locations(x_peak);
d = 0;
% Apply scaling for d value
if xor(left_intersected, right_intersected)
    if left_intersected
        d = 2 / features.distance_scaling *...
            diff([left, sensor_locations(x_peak)]);
    end
    if right_intersected
        d = 2 / features.distance_scaling *...
            diff([sensor_locations(x_peak), right]);
    end
else
    d = 1 / features.distance_scaling *...
        diff([left, right]);
end

prediction = [b, d, 0];

end

function prediction = qm_fit_location(quadrature, prediction, params)

x0 = prediction(1:2);
x_range = params.domain.x_range;
y_range = params.domain.y_range;
sensor_locations = params.sensors.locations.x;

fit_function = @(x, xdata) simulate_quadrature([x, prediction(3)], params);
options = optimoptions('lsqcurvefit',...
    'MaxIterations', params.qm.refine_iterations,...
    'FiniteDifferenceStepSize', params.qm.gradient_step,...
    'FiniteDifferenceType', 'central',...
    'StepTolerance', params.qm.step_tolerance,...
    'FunctionTolerance', params.qm.function_tolerance,...
    'OptimalityTolerance', params.qm.optimality_tolerance,...
    'Display', 'off');
prediction(1:2) = lsqcurvefit(fit_function, x0, sensor_locations, quadrature,...
                                [min(x_range), min(y_range)],... lower bound
                                [max(x_range), max(y_range)],... upper bound
                                options);

end

function prediction = qm_estimate_orientation(vx, vy, prediction, params)

wavelets = compute_wavelets_raw(prediction(1:2), params);
envelope_vx = sqrt(wavelets.even.^2 + wavelets.odd.^2);
envelope_vy = sqrt(wavelets.odd.^2 + wavelets.navelet.^2);
phi_prime_vx = unwrap(atan2(wavelets.odd, wavelets.even));
phi_prime_vy = unwrap(atan2(wavelets.navelet, wavelets.odd));

alpha_vx = real(acos(vx ./ envelope_vx));
alpha_vy = real(acos(vy ./ envelope_vy));

phi_hat_vx_plus = mod(phi_prime_vx + alpha_vx, 2*pi);
phi_hat_vx_min = mod(phi_prime_vx - alpha_vx, 2*pi);
phi_hat_vy_plus = mod(phi_prime_vy + alpha_vy, 2*pi);
phi_hat_vy_min = mod(phi_prime_vy - alpha_vy, 2*pi);

angles = [phi_hat_vx_plus phi_hat_vx_min phi_hat_vy_plus phi_hat_vy_min];
phi_hat = mod(compute_angle_median(angles), 2*pi);

prediction(3) = phi_hat;

end
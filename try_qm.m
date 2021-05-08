%% Introduction
% This shows how the quadrature predictor is used.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Get parameters
params = generate_parameters();
% params.sensors.n_sensors = 64;
% params.sensors.locations = compute_sensor_locations(params.sensors);
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
load('./config/locations/train_sources_res0.03.mat');
sources = train_sources;
clear window train_sources

%% Compute velocity
velocity = compute_sensor_velocity(sources, time, params);
velocity = apply_sensor_noise(velocity, params);
velocity = extract_signal_over_sensors(velocity, params);

%% Compute predictions
tic 
[predictions, initial_estimate] = predict_qm(velocity, params);
toc
median(compute_prediction_errors(predictions, sources, params))

% check initial orientation estimation
update_predictions = predictions;
update_predictions(:,3) = initial_estimate(:, 3);
median(compute_prediction_errors(update_predictions, sources, params))

plot_predicted_locations(predictions, sources, params);

figure
plot3(sources(:, 3), sources(:,2), predictions(:, 3), '.')
hold on
plot3(sources(:, 3), sources(:,2), initial_estimate(:, 3), '.')
grid on

%% Single source estimate
p = [0.0 0.1 0.7*pi];
velocity = compute_sensor_velocity(p, time, params);
velocity = apply_sensor_noise(velocity, params);
velocity = extract_signal_over_sensors(velocity, params);

qm = params.qm;
quadrature = sqrt(velocity.xt.^2 + 0.5 .* velocity.yt.^2);
quadrature = permute(quadrature, [1 3 2]);

[prediction, initial_estimate] = predict_qm(velocity, params);

sensor_locations = params.sensors.locations.x;

% figure(5)
% hold off
% plot(sensor_locations, quadrature, '-');
% hold on
% quadrature = simulate_quadrature(prediction, params);
% plot(sensor_locations, quadrature, '--');
% quadrature = simulate_quadrature([initial_estimate 0], params);
% plot(sensor_locations, quadrature, ':');
% 
% legend('measured', 'predicted', 'initial guess')
% xlabel('$x$ (m)')
% ylabel('weighted velocity magnitude')

% figure(6)
% hold off
% plot_predicted_locations(prediction, p, params);

figure(7)
hold off
polarplot([p(3) p(3)], [0 1]);
hold on
polarplot([prediction(3) prediction(3)], [0 1]);
polarplot([initial_estimate(3) initial_estimate(3)], [0 1]);
rticks([])
legend('actual', 'predicted', 'initial guess')

%% Compare different quadrature approaches
% 
params = generate_parameters();
% params.sensors.n_sensors = 64;
% params.sensors.locations = compute_sensor_locations(params.sensors);
load('./config/parameters/window.mat');
params.window = window;
load('./config/locations/train_sources_res0.01.mat');
sources = train_sources;
clear window train_sources

% velocity = generate_noisy_reduced_velocity(sources, params);
% save('tmp_vel.mat', 'velocity')
load('tmp_vel.mat', 'velocity')

for idx = 1:length(repeats)
    params.qm.refine_repeats = repeats(idx);
    [~, ~, ip] = predict_qm_test(velocity, params);
end

%%
repeats = 0:params.qm.refine_repeats;
median_error = zeros(1, size(ip, 1));
mean_error = zeros(1, size(ip, 1));
for rep = 1:size(ip, 1)
   errors = compute_prediction_errors(squeeze(ip(rep, :, :)), velocity.sources, params); 
   median_error(rep) = median(errors);
   mean_error(rep) = mean(errors);
   
   analyse_predictions(squeeze(ip(rep, :, :)), velocity, params);
end

table(repeats, median_error, mean_error, 'VariableNames', {'repeats', 'median', 'mean'})

figure
hold on
plot(repeats, mean_error)
plot(repeats, median_error)

%% Functions

function [predictions, estimates, ip] = predict_qm_test(velocity, params)

% Preallocation and easy access
n_samples = size(velocity.xt, 1);
velocity_xt = velocity.xt;
velocity_yt = velocity.yt;
ip  = zeros(params.qm.refine_repeats+1, n_samples, 3);

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
    sample_ip = ip(:, sample_idx, :);
    
    sample_ip(1, :) = qm_feature_estimate(sample_quadrature_norm, features, params.sensors.locations.x);
    sample_ip(1, :) = qm_estimate_orientation(sample_vx, sample_vy, sample_ip(1, :), params);

    for repeat_idx = (1:params.qm.refine_repeats)+1
       sample_ip(repeat_idx, :) = qm_fit_location(sample_quadrature, sample_ip(repeat_idx-1, :), params);
       sample_ip(repeat_idx, :) = qm_estimate_orientation(sample_vx, sample_vy, sample_ip(repeat_idx, :), params);
    end
    
    ip(:, sample_idx, :) = sample_ip;
end

estimates = ip(1, :, :);
predictions = ip(end, :, :);

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


function [mean_error, median_error] = analyse_predictions(predictions, velocity, params)

cell_size=0.02;

location_errors = compute_location_error(predictions, velocity.sources);
orientation_errors = compute_orientation_error(predictions, velocity.sources);
prediction_error = compute_prediction_errors(predictions, velocity.sources, params);

% Compute spatial bins
bins = compute_bins(velocity.sources, cell_size);

predictions = table(repmat('mlp', size(predictions, 1), 1),...
                    bins, predictions, location_errors,...
                    orientation_errors, prediction_error,...
                    'VariableNames', {'method', 'source_bin', 'prediction',...
                                      'location_error', 'orientation_error', 'error'});
[xy_bins, x, y] = findgroups(...
    predictions.source_bin(:,1),...
    predictions.source_bin(:,2));
median_location_error = splitapply(@median, predictions.location_error, xy_bins);
median_orientation_error = splitapply(@median, abs(predictions.orientation_error), xy_bins);
median_prediction_error = splitapply(@median, predictions.error, xy_bins);

spatial_error = table(x, y, median_location_error, median_orientation_error, median_prediction_error);
x = unique(x);
y = unique(y);
[X, Y] = meshgrid(x, y);

conf.thresholds = [0 0.01, 0.03, 0.05, 0.09];
conf.cell_size = cell_size;

figure
hold off
data = reshape(spatial_error.median_location_error, [length(y), length(x)]);
contourf(X, Y, data, conf.thresholds(1:end))
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size/2 max(y)])
xlabel('$x$ (m)')
ylabel('$y$ (m)')
colormap(flip(winter(1024)))
caxis([0 conf.thresholds(end)])
daspect([1 1 1])
drawnow
pause(0.2)
post_process_figure(1, 0.5, [1 1], [0 0])
set(gca, 'tickdir', 'out')


figure
hold off
data = reshape(spatial_error.median_orientation_error / pi, [length(y), length(x)]);
contourf(X, Y, data, conf.thresholds(1:end))
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size/2 max(y)])
xlabel('$x$ (m)')
ylabel('$y$ (m)')
colormap(flip(autumn(1024)))
caxis([0 conf.thresholds(end)])
daspect([1 1 1])
drawnow
pause(0.2)
post_process_figure(1, 0.5, [1 1], [0 0])
set(gca, 'tickdir', 'out')

mean_error = mean(prediction_error);
median_error = median(prediction_error);
disp(' ')
disp(['Median Absolute Validation Error: ', num2str(median_error)]);
disp(['Mean Absolute Validation Error: ', num2str(mean_error)]);
disp(' ')

end

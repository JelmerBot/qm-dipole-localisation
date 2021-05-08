%% Introduction
% This shows how the cwt predictor is used.

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
load('./config/locations/train_sources_res0.09.mat');
validate_sources = train_sources;
load('./config/locations/train_sources_res0.03.mat');
clear window

%% Compute train velocity
% CWT needs even, odd and navelet at the required distances
wavelets = compute_wavelets(train_sources, params);

validate_velocity = compute_sensor_velocity(validate_sources, time, params);
validate_velocity = apply_sensor_noise(validate_velocity, params);
validate_velocity = extract_signal_over_sensors(validate_velocity, params);

%% Compute predictions (no normalisation)
% maybe normalise components?
tic 
predictions = predict_cwt(validate_velocity, wavelets, params);
toc
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%% Compute predictions (component normalisation)
% maybe normalise components?
velocity = normalise_components(validate_velocity, params);
tic 
predictions = predict_cwt(velocity, wavelets, params);
toc
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%%

filter = validate_sources(:, 1) > -0.2 & validate_sources(:,1) < 0.2 &...
         validate_sources(:, 2) < 0.3;

dif_loc = predictions(filter, 1:2) - validate_sources(filter, 1:2);
dif_or = compute_angle_difference(validate_sources(filter, 3), predictions(filter, 3)) / pi;

figure
subplot(1,3,1)
plot(validate_sources(filter,3) / pi, dif_loc(:,1), '.')
xlabel('Orientation source')
ylabel('X difference')
subplot(1,3,2)
plot(validate_sources(filter,3) / pi, dif_loc(:,2), '.')
xlabel('Orientation source')
ylabel('Y difference')
subplot(1,3,3)
plot(validate_sources(filter,3) / pi, dif_or(:,1), '.')
xlabel('Orientation source')
ylabel('Phi difference')

%% Compute one source prediction
% params.sensors.input_mode = 'x+y';
p = [0.1 0.2 0.2*pi];
velocity = compute_sensor_velocity(p, 0, params);

% Convinience variables
sensors = params.sensors;
n_samples = size(velocity.xt, 1);
n_wavelets = size(wavelets.even, 1);

% Apply input mode
[vx, vy] = cwt_input_mode_velocity(velocity, params);
vx = permute(vx, [1 3 2]);
vy = permute(vy, [1 3 2]);
even = wavelets.even;
odd = wavelets.odd;
navelet = wavelets.navelet;

% Compute cwt coefficients
scaling = 1 ./ sqrt(abs(sensors.y_value - wavelets.sources(:,2)))';
x_even = zeros(n_samples, n_wavelets);
x_odd = zeros(n_samples, n_wavelets);
y_odd = zeros(n_samples, n_wavelets);
y_nav = zeros(n_samples, n_wavelets);
% vx is measured by the sensors?
if ~strcmp(params.sensors.input_mode, 'y')
    x_even = scaling .* (vx * even');
    x_odd =  scaling .* (vx * odd'); 
end
% vy is measured by the sensors?
if ~strcmp(params.sensors.input_mode, 'x')
    y_nav =  scaling .* (vy * navelet');
    y_odd =  scaling .* (vy * odd');
end
magnitude = sqrt(x_even.^2 + y_nav.^2 + x_odd.^2 + y_odd.^2);
Z = magnitude ./ max(magnitude(:));
X = wavelets.sources(:,1);
Y = wavelets.sources(:,2);

% Fit gaussian model to magnitude
thresholds = [params.cwt.threshold_min, params.cwt.threshold_max];
filter = Z > thresholds(1) & Z < thresholds(2);
n_sources = sum(filter(:));
if n_sources <= 0
    warning('CWT:filterExhausted', ...
        ['No coefficient values passed the threshold filter. '...
         'Used the maximum value instead.'])
    [~, i] = max(X(:));
    p_hat = [X(i), Y(i)];
elseif n_sources < 7 % The number of variables that is fitted
    warning('CWT:filterTooRestrictive', ...
        ['Not enough points passed the threshold filter to fit a Gaussian.\n',...
         'Used the maximum value within the filter instead.'])
    [~, i] = max(Z(filter));
    X_filter = X(filter);
    Y_filter = Y(filter);
    p_hat = [X_filter(i), Y_filter(i)];
else
    res = fmgaussfit(X(filter), Y(filter), Z(filter)', params);
    p_hat = res(5:6);
end

% Compute orientation
p_hat_wavelets = compute_wavelets(p_hat, params);
p_hat_scaling = 1 ./ sqrt(abs(sensors.y_value - p_hat_wavelets.sources(:,2)));
% vx is measured by the sensors?
if ~strcmp(params.sensors.input_mode, 'y')
    p_hat_x_even = p_hat_scaling .* sum(vx .* p_hat_wavelets.even, 2);
    p_hat_x_odd =  p_hat_scaling .* sum(vx .* p_hat_wavelets.odd, 2);
    phi_x = mod(atan2(params.cwt.c_x * p_hat_x_odd, p_hat_x_even), 2*pi); 
end
% vy is measured by the sensors?
if ~strcmp(params.sensors.input_mode, 'x')
    p_hat_y_nav =  p_hat_scaling .* sum(vy .* p_hat_wavelets.navelet, 2);
    p_hat_y_odd =  p_hat_scaling .* sum(vy .* p_hat_wavelets.odd, 2);
    phi_y = mod(atan2(params.cwt.c_y * p_hat_y_nav, p_hat_y_odd), 2*pi);
end

switch params.sensors.input_mode
    case 'x'
        phi = phi_x;
    case 'y'
        phi = phi_y;
    otherwise % x+y and x|y
        phi = compute_angle_average([phi_x, phi_y]);
end
p_hat = [p_hat phi];

% Plot values
figure(1)
plot_cwt(X, Y, Z, x_even, y_nav, x_odd, y_odd, p, p_hat, params)

%% Plot the fitted gaussian

figure(2)
hold off
% Coefficients
a1 = plot3(X(filter), Y(filter), Z(filter), '.');
hold on
a2 = plot3(X(~filter), Y(~filter), Z(~filter), '.');
% Source location
a3 = plot3(p(1), p(2), 0, 'kx');
plot3(p(1), p(2), 1, 'kx')
plot3([p(1) p(1)], [p(2) p(2)], [0, 1], 'k-')
% Fitted Gaussian
[X_fit, Y_fit] = meshgrid(-0.5:0.02:0.5, 0:0.02:0.5);
xy = {X_fit(:),Y_fit(:)};
fitted_gaus = gaussian2D(res, xy);
a5 = surf(X_fit, Y_fit, reshape(fitted_gaus, size(X_fit)));
alpha 0.75
a4 = plot3(p_hat(1), p_hat(2), 0, 'go');
plot3(p_hat(1), p_hat(2), 1, 'go')
plot3([p_hat(1) p_hat(1)], [p_hat(2) p_hat(2)], [0, 1], 'g-')
% Layout
grid on
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$Wv$ (normalised)')
legend([a1 a2 a3 a4 a5], 'used Coefficients', 'unused Coefficients',...
             'source location', 'estimated location', 'fitted Gaussian')
view(-39.9, 10.8)
print('cwt_fitting.jpg', '-djpeg', '-r600')

%% Create location and orientation error plots 
% checks validation of the parameters
params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
load('./config/locations/test_sources.mat', 'test_sources')
load('./config/locations/train_sources_res0.09.mat', 'train_sources')
load('config/parameters/cwt_res0.09.mat', 'cwt')
params.cwt = cwt;
cell_size = 0.02;

wavelets = compute_wavelets(train_sources, params);
% velocity = generate_noisy_reduced_velocity(test_sources, params);
% save('tmp_vel.mat', 'velocity');
load('tmp_vel.mat', 'velocity');

predictions = predict_cwt(velocity, wavelets, params);

location_errors = compute_location_error(predictions, velocity.sources);
orientation_errors = compute_orientation_error(predictions, velocity.sources);
prediction_error = compute_prediction_errors(predictions, velocity.sources, params);

% Compute spatial bins
bins = compute_bins(velocity.sources, cell_size);

predictions = table(repmat('cwt', size(predictions, 1), 1),...
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
%%
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
% post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
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
% post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
set(gca, 'tickdir', 'out')

figure
hold off
data = reshape(spatial_error.median_prediction_error, [length(y), length(x)]);
contourf(X, Y, data, conf.thresholds(1:end))
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size/2 max(y)])
xlabel('$x$ (m)')
ylabel('$y$ (m)')
colormap(flip(hot(1024)))
caxis([0 conf.thresholds(end)])
daspect([1 1 1])
% post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
set(gca, 'tickdir', 'out')

figure
plot_predicted_locations(predictions.prediction(1:100,:), velocity.sources(1:100,:), params)

%% Introduction
% This shows how the gn predictor is used.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();
gcp;

%% Get parameters
params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
load('./config/locations/train_sources_res0.09.mat');
sources = train_sources;
clear window train_sources

%% Compute velocity
velocity = compute_sensor_velocity(sources, time, params);
velocity = apply_sensor_noise(velocity, params);
velocity = extract_signal_over_sensors(velocity, params);
velocity = apply_input_mode(velocity, params);

%% Call the predictor
tic
predictions = predict_gn(velocity, params);
toc
median(compute_prediction_errors(predictions, sources, params))
plot_predicted_locations(predictions, sources, params);

%% Plot single estimate

[vel, idx] = select_velocity(velocity, [0, 1], []);
[p, estimates] = predict_gn(vel, params);

lim = find(isnan(estimates(1, :, 1)), 1) - 1;
if isempty(lim)
    lim = params.gn.max_iterations;
end
figure(2)
hold off
plot3(estimates(1, 1:lim, 1), estimates(1, 1:lim, 2), estimates(1, 1:lim, 3) / pi, '.')
hold on
plot3(estimates(1, lim, 1), estimates(1, lim, 2), estimates(1, lim, 3) / pi, 'go')
plot3(estimates(1, 1, 1), estimates(1, 1, 2), estimates(1, 1, 3) / pi, 'ro')
plot3(sources(idx,1), sources(idx, 2), sources(idx, 3) / pi, 'kx')
grid on
xlim(params.domain.x_range)
ylim(params.domain.y_range)
zlim(params.domain.azimuth_range / pi)
view(3)

figure(3)
colors = get(gca, 'colororder');
% Noiseless expected velocity
tmp = compute_sensor_velocity(vel.sources, 0, params);
hold off
plot(squeeze(abs(tmp.xt)), ':', 'color', colors(1, :));
hold on
plot(squeeze(abs(tmp.yt)), ':', 'color', colors(2, :));
% Noisy reconstructed velocity
plot(squeeze(abs(vel.xt)), '-', 'color', colors(1, :));
plot(squeeze(abs(vel.yt)), '-', 'color', colors(2, :));
% Solution by GN
tmp = compute_sensor_velocity(estimates(1, lim, :), 0, params); 
plot(squeeze(abs(tmp.xt)), '--', 'color', colors(1, :));
plot(squeeze(abs(tmp.yt)), '--', 'color', colors(2, :));

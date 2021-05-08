%% Introduction
% This shows how the lcmv predictor is used.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Get parameters
params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
load('./config/locations/train_sources_res0.09.mat');
validate_sources = train_sources;
load('./config/locations/train_sources_res0.03.mat');
clear window

%% Compute train velocity
train_velocity = compute_sensor_velocity(train_sources, time, params);
train_velocity = apply_sensor_noise(train_velocity, params);
train_velocity = extract_signal_over_sensors(train_velocity, params);
validate_velocity = compute_sensor_velocity(validate_sources, time, params);
validate_velocity = apply_sensor_noise(validate_velocity, params);

%% Call the predictor no normalisation

train = apply_input_mode(train_velocity, params);
validate = apply_input_mode(validate_velocity, params);

predictions = predict_lcmv(validate, train, params);
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%% Call the predictor input normalisation

train = apply_input_mode(train_velocity, params);
train = normalise_input(train);
validate = apply_input_mode(validate_velocity, params);
% validate = normalise_input(validate);

predictions = predict_lcmv(validate, train, params);
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%% Call the predictor component normalisation

train = normalise_components(train_velocity, params);
train = apply_input_mode(train, params);
% validate = normalise_components(validate_velocity, params);
validate = apply_input_mode(validate, params);

predictions = predict_lcmv(validate, train, params);
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%% Plot single estimate

p = [0.05 0.1 1.7*pi];
vel = compute_sensor_velocity(p, time, params);
vel = apply_sensor_noise(vel, params);
show_velocity = extract_signal_over_sensors(vel, params);
train = train_velocity;

n_modes = length(params.experiment.input_modes);
cnt = 1;

for mode = params.experiment.input_modes
    params.sensors.input_mode = mode{1};
    
    v = apply_input_mode(vel, params);
    t = apply_input_mode(train, params);
    t = normalise_input(t);
    
    [p1, E1] = predict_lcmv(v, t, params);
    [v1, i1] = max(E1(:));
    [vk, i] = maxk(E1(:), 20);
    
    figure(4)
    subplot(2, n_modes, cnt)
    hold off
    scatter3(train_sources(i, 1),...
            train_sources(i, 2),...
            train_sources(i, 3) / pi,...
            5,...
            log(vk ./ v1), 'filled');
    hold on
    plot3(p1(1), p1(2), p1(3)/pi, 'go')
    plot3(v.sources(1), v.sources(2), v.sources(3)/pi, 'kx')
    zlim([0, 2])
    grid on
    view(3)
    
    subplot(2, n_modes, n_modes + cnt)
    colors = get(gca, 'colororder');
    % Noisy actual location
    hold off
    plot(squeeze(show_velocity.xt), '-', 'color', colors(1, :));
    hold on
    plot(squeeze(show_velocity.yt), '-', 'color', colors(2, :));
    % Noiseless actual location
    tmp = compute_sensor_velocity(show_velocity.sources, 0, params);
    plot(squeeze(tmp.xt), ':', 'color', colors(1, :));
    plot(squeeze(tmp.yt), ':', 'color', colors(2, :));
    
    % Noiseless predicted location
    tmp = compute_sensor_velocity(p1, 0, params);
    plot(squeeze(tmp.xt), '--', 'color', colors(1, :));
    plot(squeeze(tmp.yt), '--', 'color', colors(2, :));
    
    cnt = cnt + 1;
end
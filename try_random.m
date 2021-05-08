%% Introduction
% This shows how the random predictor is used.

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
load('./config/locations/train_sources_res0.06.mat');
sources = train_sources;
clear window train_sources

%% Compute velocity
% Only really need the size...
velocity = compute_sensor_velocity(sources, time, params);
velocity = apply_sensor_noise(velocity, params);
velocity = extract_signal_over_sensors(velocity, params);
velocity = apply_input_mode(velocity, params);

%% Call the predictor
predictions = predict_random(velocity, params);
mean(compute_prediction_errors(predictions, sources, params))
plot_predicted_locations(predictions, sources, params);

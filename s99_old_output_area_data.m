%% Introduction
% Demonstrates the effect of the area in which source states are located on
% the error distrubtion of GN.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();
gcp;

%% Load parameters

params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
load('./config/parameters/gn_input_mode_1.mat', 'gn');
params.gn = gn;

%% Evaluate small area

load('./config/locations/train_sources_res0.03.mat', 'train_sources');
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);
predictions = predict_gn(velocity, params);

small_errors = compute_prediction_errors(predictions, train_sources, params);

%% Evaluate larger area

% Update parameters
params.domain.x_range = [-1 1];
params.domain.y_range = [0 1] + params.source.radius;
params.domain.min_source_distance = ...
        params.experiment.train_source_distances(2);
    
% Generate sources
train_sources = compute_sampling_locations(params);
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);
predictions = predict_gn(velocity, params);

large_errors = compute_prediction_errors(predictions, train_sources, params);

%% Create plot

figure()
subplot(121)
histogram(small_errors)
subplot(122)
histogram(large_errors)

%% Save data

res = table([small_errors; large_errors],...
            [repmat('small', size(small_errors)); 
             repmat('large', size(large_errors))]);
res.Properties.VariableNames = {'error', 'area'};

writetable(res, 'area.csv')
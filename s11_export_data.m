%% Introduction
% Combines all predictions into two csv files.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

% Processing variables
cell_size = 0.02;

% Make sure to update these when you've ran the simulation!
input_folder_data = './data/20_05_21/';
input_folder_sensors = './data/20_05_21/';
output_folder = './data/final_csv/';

%% Merge output test_methods_data
files = dir([input_folder_data, 'predictions_res*.mat']);
load([input_folder_data, 'params_res0.01.mat'], 'params');

merged = [];
for idx = 1:length(files)
    load([input_folder_data, files(idx).name])
    name = files(idx).name;
    pos = strfind(name, '_res');
    res = str2double(name(pos+4:end-4));
    predictions.resolution = repmat(res, height(predictions), 1);
    merged = [merged; predictions];
    clear predictions
end
predictions = merged;
predictions.Properties.VariableNames{1} = 'method';
predictions.method = categorical(predictions.method);

% Add location bins
predictions.source_bin = compute_bins(predictions.source, cell_size);

% Add errors
predictions.prediction_error = compute_prediction_errors(predictions.prediction, predictions.source, params);
predictions.location_error = compute_location_error(predictions.prediction, predictions.source);
predictions.orientation_error = compute_orientation_error(predictions.prediction, predictions.source);
predictions.orientation_error_no_phase = -abs(predictions.orientation_error - pi/2) + pi/2;

% Split components for csv write
predictions.x_predicted = predictions.prediction(:, 1);
predictions.y_predicted = predictions.prediction(:, 2);
predictions.orientation_predicted = predictions.prediction(:, 3);

predictions.x_source = predictions.source(:, 1);
predictions.y_source = predictions.source(:, 2);
predictions.orientation_source = predictions.source(:, 3);

predictions.x_bin = predictions.source_bin(:, 1);
predictions.y_bin = predictions.source_bin(:, 2);
predictions.orientation_bin = predictions.source_bin(:, 3);

% drop vector columns
predictions = removevars(predictions, {'prediction', 'source', 'source_bin'});

% write as csv
writetable(predictions, [output_folder, 'predictions_resolution.csv'])

%% Merge output test_methods_input

files = dir([input_folder_sensors, 'predictions_input_mode*.mat']);
load([input_folder_sensors, 'params_input_mode_1.mat'], 'params');

merged = [];
for idx = 1:length(files)
    load([input_folder_sensors, files(idx).name])
    name = files(idx).name;
    mode = round(str2double(name(end-4:end-4)));
    predictions.input_mode = repmat(params.experiment.input_modes(mode), height(predictions), 1);
    merged = [merged; predictions];
    clear predictions
end
predictions = merged;
predictions.Properties.VariableNames{1} = 'method';
predictions.input_mode = categorical(predictions.input_mode);
predictions.method = categorical(predictions.method);

% Fix names
% predictions.method = upper(arrayfun(@char, predictions.method, 'UniformOutput', false));
% predictions.method = categorical(predictions.method);
% u_method = unique(predictions.method);
% n_method = u_method;
% n_method(n_method == 'RAND') = 'RND';
% predictions.method = renamecats(predictions.method, arrayfun(@char, u_method, 'UniformOutput', false), arrayfun(@char, n_method, 'UniformOutput', false));

% Add location bins
predictions.source_bin = compute_bins(predictions.source, cell_size);

% Add errors
predictions.prediction_error = compute_prediction_errors(predictions.prediction, predictions.source, params);
predictions.location_error = compute_location_error(predictions.prediction, predictions.source);
predictions.orientation_error = compute_orientation_error(predictions.prediction, predictions.source);
predictions.orientation_error_no_phase = -abs(predictions.orientation_error - pi/2) + pi/2;

% Split components for csv write
predictions.x_predicted = predictions.prediction(:, 1);
predictions.y_predicted = predictions.prediction(:, 2);
predictions.orientation_predicted = predictions.prediction(:, 3);

predictions.x_source = predictions.source(:, 1);
predictions.y_source = predictions.source(:, 2);
predictions.orientation_source = predictions.source(:, 3);

predictions.x_bin = predictions.source_bin(:, 1);
predictions.y_bin = predictions.source_bin(:, 2);
predictions.orientation_bin = predictions.source_bin(:, 3);

% drop vector columns
predictions = removevars(predictions, {'prediction', 'source', 'source_bin'});

% write as csv
writetable(predictions, [output_folder, 'predictions_input_mode.csv'])


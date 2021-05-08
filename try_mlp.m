%% Introduction
% This shows how the mlp predictor is used.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% params
params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
time = compute_time(params);
clear window
load('./config/parameters/elm_input_mode_1.mat', 'elm');
load('./config/parameters/mlp_input_mode_1.mat', 'mlp');
params.mlp = mlp;
params.elm = elm;
params.mlp.execution_environment = 'gpu';
params.mlp.verbose = true;

%% Determine sizes of layers
elm_total_weights = 16 * elm.n_nodes +  elm.n_nodes * 3
elm_learned_weights = elm.n_nodes * 3
mlp_total_weights = sum([16 mlp.layer_sizes(1:end-1)] .* [mlp.layer_sizes(2:end) 3])

ratio_learned = mlp_total_weights / elm_learned_weights
ratio_total= mlp_total_weights / elm_total_weights

%% Compute velocity

res = 0.01;
load(['./config/locations/train_sources_res', num2str(res) '.mat'], 'train_sources');
% velocity = generate_noisy_reduced_velocity(train_sources, params);
% save('tmp_vel.mat', 'velocity');
load('tmp_vel.mat', 'velocity');
velocity = apply_input_mode(velocity, params);

n_fold = params.experiment.n_fold;
n_samples = size(velocity.input, 1);
partitions = crossvalind('Kfold', n_samples, n_fold);
validation_velocity = select_velocity_indices(velocity, partitions == 1);
train_velocity = select_velocity_indices(velocity, partitions ~= 1);

load('./config/locations/test_sources.mat', 'test_sources');
% test_velocity = generate_noisy_reduced_velocity(test_sources, params);
% save('tmp_vel_t.mat', 'test_velocity');
load('tmp_vel_t.mat', 'test_velocity');
test_velocity = apply_input_mode(test_velocity, params);

%% Evaluate MLP

tic
mlp = train_mlp(validation_velocity, train_velocity, params);
disp(['Training in: ', duration2str(toc)])
save('tmp_mlp.mat', 'mlp')
% load('tmp_mlp.mat', 'mlp')

tic;
predictions = predict_mlp(test_velocity, mlp, params);
disp(['Predictions in: ', duration2str(toc)])

analyse_predictions(predictions, test_velocity,  params)

%% MSE performance ReLU (no dropout no batchnorm)

params.mlp.error_metric = 'mean_squared_error';
params.mlp.use_dropout = false;
params.mlp.use_batch_normalisation = false;

validate_idx = partitions == 1;
train_idx = ~validate_idx;

train = select_velocity_indices(velocity, train_idx);
validate = select_velocity_indices(velocity, validate_idx);

[predictions, ~] = train_and_predict_mlp(validate, train, params);
analyse_predictions(predictions, validate, params);

%% MAE performance ReLU (no dropout no batchnorm)

params.mlp.error_metric = 'mean_absolute_error';
params.mlp.use_dropout = false;
params.mlp.use_batch_normalisation = false;

validate_idx = partitions == 1;
train_idx = ~validate_idx;

train = select_velocity_indices(velocity, train_idx);
validate = select_velocity_indices(velocity, validate_idx);

[predictions, ~] = train_and_predict_mlp(validate, train, params);
analyse_predictions(predictions, validate, params);

%% MedianAE performance ReLU (no dropout no batchnorm)

params.mlp.error_metric = 'median_absolute_error';
params.mlp.use_dropout = false;
params.mlp.use_batch_normalisation = false;

validate_idx = partitions == 1;
train_idx = ~validate_idx;

train = select_velocity_indices(velocity, train_idx);
validate = select_velocity_indices(velocity, validate_idx);

[predictions, ~] = train_and_predict_mlp(validate, train, params);
analyse_predictions(predictions, validate, params);

%% MAE performance ReLU (no dropout no batchnorm) no pretraining

params.mlp.error_metric = 'mean_absolute_error';
params.mlp.use_dropout = false;
params.mlp.use_batch_normalisation = false;
params.mlp.training_levels = Inf;

validate_idx = partitions == 1;
train_idx = ~validate_idx;

train = select_velocity_indices(velocity, train_idx);
validate = select_velocity_indices(velocity, validate_idx);

[predictions, ~] = train_and_predict_mlp(validate, train, params);
analyse_predictions(predictions, validate, params);

%% Helpers

function analyse_predictions(predictions, velocity, params)

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
post_process_figure(1, 0.5, [1 1], [0 0])
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
post_process_figure(1, 0.5, [1 1], [0 0])
set(gca, 'tickdir', 'out')

figure
plot_predicted_locations(predictions.prediction(1:100,:), velocity.sources(1:100,:), params)

disp(' ')
disp(['Median Absolute Validation Error: ', num2str(median(prediction_error))]);
disp(['Mean Absolute Validation Error: ', num2str(mean(prediction_error))]);
disp(' ')

end
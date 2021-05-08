%% Introduction
% This shows how the elm predictor is used.

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
load('./config/locations/train_sources_res0.01.mat');
load('./config/locations/test_sources.mat');
clear window
load('./config/parameters/elm_input_mode_4.mat', 'elm');
params.elm = elm;
clear elm

params.sensors.input_mode = 'x|y';

%% Compute train velocity

train_velocity = generate_noisy_reduced_velocity(train_sources, params);
train_velocity = apply_input_mode(train_velocity, params);
train_velocity = normalise_input(train_velocity);

%% Compute validate velocity

validate_velocity = generate_noisy_reduced_velocity(test_sources, params);
validate_velocity = apply_input_mode(validate_velocity, params);
validate_velocity = normalise_input(validate_velocity);

%% Evaluate
tic
elm = train_elm(normalise_input(train_velocity), params);
disp(duration2str(toc))
tic
predictions = predict_elm(normalise_input(validate_velocity), elm, params);
disp(duration2str(toc))

validate_error = median(compute_prediction_errors(predictions, test_sources, params))
sum(isnan(predictions))

%% Call the predictor with tanh

errors = zeros(1, 100);
for idx = 1:100
params.elm.f_internal = @tanh;
[predictions, training_error, train_predictions] = train_and_predict_elm(validate_velocity, train_velocity, params);
validate_error = median(compute_prediction_errors(predictions, validate_sources, params));
errors(idx) = validate_error;

end

figure
plot(errors);
title('Tanh')
drawnow
% disp(['validation error: ', num2str(validate_error)]);
% disp(['traiing error: ', num2str(training_error)]);
figure
plot_predicted_locations(predictions, validate_sources, params);
figure
plot_predicted_locations(train_predictions, train_sources, params);

%% Call the predictor with relu

errors = zeros(1, 100);
for idx = 1:100
params.elm.f_internal = @relu;
[predictions, training_error, train_predictions] = train_and_predict_elm(validate_velocity, train_velocity, params);
validate_error = median(compute_prediction_errors(predictions, validate_sources, params));
errors(idx) = validate_error;

end
figure
plot(errors);
title('ReLu')
drawnow
% disp(['validation error: ', num2str(validate_error)]);
% disp(['training error: ', num2str(training_error)]);
figure
plot_predicted_locations(predictions, validate_sources, params);
figure
plot_predicted_locations(train_predictions, train_sources, params);

%% Call the predictor with dCaAp

errors = zeros(1, 100);
for idx = 1:100
params.elm.f_internal = @dCaAp;
[predictions, training_error, train_predictions] = train_and_predict_elm(validate_velocity, train_velocity, params);
validate_error = median(compute_prediction_errors(predictions, validate_sources, params));
errors(idx) = validate_error;

end

figure
plot(errors);
title('dCaAp')
drawnow
% disp(['validation error: ', num2str(validate_error)]);
% disp(['training error: ', num2str(training_error)]);
figure
plot_predicted_locations(predictions, validate_sources, params);
figure
plot_predicted_locations(train_predictions, train_sources, params);

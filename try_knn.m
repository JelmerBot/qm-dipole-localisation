%% Introduction
% This shows how the knn predictor is used.

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
load('./config/locations/train_sources_res0.03.mat');
validate_sources = train_sources;
load('./config/locations/train_sources_res0.05.mat');
clear window

%% Compute train velocity
train_velocity = compute_sensor_velocity(train_sources, time, params);
train_velocity = apply_sensor_noise(train_velocity, params);
train_velocity = extract_signal_over_sensors(train_velocity, params);
train_velocity = apply_input_mode(train_velocity, params);
train_velocity = normalise_input(train_velocity);

%% Compute validate velocity
validate_velocity = compute_sensor_velocity(validate_sources, time, params);
validate_velocity = apply_sensor_noise(validate_velocity, params);
validate_velocity = extract_signal_over_sensors(validate_velocity, params);
validate_velocity = apply_input_mode(validate_velocity, params);
validate_velocity = normalise_input(validate_velocity);

%% Call the predictor
predictions = predict_knn(validate_velocity, train_velocity, params);
median(compute_prediction_errors(predictions, validate_sources, params))
plot_predicted_locations(predictions, validate_sources, params);

%% Plot single estimate

[vel, idx] = select_velocity(validate_velocity, [0.1, 0.3], [1.5*pi]);
[p, indices] = predict_knn(vel, train_velocity, params);

figure(2)
hold off
plot3(train_velocity.sources(:,1),...
      train_velocity.sources(:,2),...
      train_velocity.sources(:,3)/pi, '.')
hold on
plot3(vel.sources(1), vel.sources(2), vel.sources(3)/pi, 'kx')
plot3(p(1), p(2), p(3)/pi, 'go')
plot3(train_velocity.sources(indices, 1),...
     train_velocity.sources(indices, 2),...
     train_velocity.sources(indices,3)/pi, 'r.', 'MarkerSize', 10)

x_range = params.domain.x_range;
y_range = params.domain.y_range;

axis image
xlim([-0.55 0.55])
ylim([0, 0.5])
zlim([0, 2])
grid on

xlabel('$x$~(m)')
ylabel('$y$~(m)')
zlabel('$\varphi$~(rad)')

drawnow
post_process_figure(1, 0.87, [1.05 1], [0.2 0.2])
curunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
cursize = get(gca, 'Position');
set(gca, 'Units', curunits);
pt_sz = diff(xlim) / cursize(3);

figure(3)
colors = get(gca, 'colororder');
% Noiseless expected velocity
tmp = compute_sensor_velocity(vel.sources, 0, params);
hold off
plot(squeeze(tmp.xt), '-.', 'color', colors(1, :));
hold on
plot(squeeze(tmp.yt), '-.', 'color', colors(2, :));
% Noisy reconstructed velocity
plot(squeeze(vel.xt), '-', 'color', colors(1, :));
plot(squeeze(vel.yt), '-', 'color', colors(2, :));
% Predicted
tmp = compute_sensor_velocity(p, 0, params);
plot(squeeze(tmp.xt), '--', 'color', colors(1, :));
plot(squeeze(tmp.yt), '--', 'color', colors(2, :));
% Plot k neighbors
plot(squeeze(train_velocity.xt(indices, :, :))', ':', 'color', colors(1,:))
plot(squeeze(train_velocity.yt(indices, :, :))', ':', 'color', colors(2,:))
xlabel('sensor')
ylabel('normalised velocity')
legend('$v_x$', '$v_y$', 'Location', 'northwest')
post_process_figure(1, 0.6, [1.05 1], [0.2 0.4])

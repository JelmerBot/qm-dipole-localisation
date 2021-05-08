%% Introduction
% Demonstrates the difference between window functions and the way time is 
% sampled used in the FFT when averaging out the time-dimension from the 
% velocity measurements.

% This file was used to check our chosen approaches and outputs hamming
% window regardless of the reconstruction errors encountered.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Get parameters
params = generate_parameters();
params.sensors.n_sensors = 8;
params.sensors.locations = compute_sensor_locations(params.sensors);
params.output_folder = './config/parameters';
if ~exist(params.output_folder, 'dir')
    mkdir(params.output_folder)
end

time = compute_time(params);

load('./config/locations/train_sources_res0.03.mat');
sources = train_sources;
clear train_sources

% Prototype wavelet
control = compute_sensor_velocity(sources, 0, params);

% Noisy velocity
noiseless = compute_sensor_velocity(sources, time, params);
noisy = apply_sensor_noise(noiseless, params);

%% Reconstructions
windows = {@rectwin, @hann, @hamming, @flattopwin};
n_windows = length(windows);
errors = zeros(1, n_windows);
raw_errors = zeros(size(control.xt, 1), n_windows);

[~, idx] = select_velocity(control, [0.05, 0.4], []);

figure(1)
colors = get(gca, 'colororder');
cnt = 1;
for window = windows
    
    params.window = window{1};
    noisy_period = extract_signal_over_sensors(noisy, params);    
    
    subplot(1, n_windows, cnt)
    % Control
    hold off
    plot(squeeze(control.xt(idx, :, :)), '-', 'color', colors(1, :));
    hold on
    plot(squeeze(control.yt(idx, :, :)), '-', 'color', colors(2, :));
    % Noiseless
    plot(squeeze(noisy_period.xt(idx, :, :)), '--', 'color', colors(1, :));
    plot(squeeze(noisy_period.yt(idx, :, :)), '--', 'color', colors(2, :));
    raw_errors(:, cnt) = sum(...
        (permute(cat(3, noisy_period.xt, noisy_period.yt), [1 3 2]) -...
         permute(cat(3, control.xt, control.yt), [1 3 2])).^2, 2);
    errors(cnt) = mean(raw_errors(:, cnt));
    cnt = cnt + 1;
end

figure
for idx = 1:n_windows
   subplot(n_windows, 1, idx) 
   histogram(raw_errors(:, idx), 500)
   xlim([0, 1e-8])
   title(func2str(windows{idx}))
end

disp('Window variants: ')
arrayfun(@(x) disp([func2str(windows{x}) ': ' num2str(mean(raw_errors(:, x))), ' m/s']), 1:n_windows);

[~, i] = min(errors);
% window = windows{i};
window = windows{3}; % Use hamming, even if Rect works better
save(fullfile(params.output_folder, 'window.mat'), 'window')

%% Compare time methods
sampling_rate = params.sensors.sampling_rate;
duration = params.domain.measurement_duration;
period = 1 / sampling_rate;
total_length = sampling_rate * params.domain.measurement_duration;
times = {...
    @() (1:total_length) * period,...
    @() (0:total_length-1) * period,...
    @() linspace(0, duration, total_length),...
    @() linspace(0, duration, total_length + 1),...
};
n_times = length(times);
errors = zeros(1, n_times);

[~, idx] = select_velocity(control, [0.05, 0.4], []);

figure(2)
colors = get(gca, 'colororder');
cnt = 1;
for time_fun = times
    
    time = time_fun{1}();
    noiseless = compute_sensor_velocity(sources, time, params);
    rng(0); noisy = apply_sensor_noise(noiseless, params);
    noisy_period = extract_signal_over_sensors(noisy, params);    
    
    subplot(1, n_times, cnt)
    % Control
    hold off
    plot(squeeze(control.xt(idx, :, :)), '-', 'color', colors(1, :));
    hold on
    plot(squeeze(control.yt(idx, :, :)), '-', 'color', colors(2, :));
    % Noiseless
    plot(squeeze(noisy_period.xt(idx, :, :)), '--', 'color', colors(1, :));
    plot(squeeze(noisy_period.yt(idx, :, :)), '--', 'color', colors(2, :));
    errors(cnt) = mean(sum(...
        (permute(cat(3, noisy_period.xt, noisy_period.yt), [1 3 2]) -...
         permute(cat(3, control.xt, control.yt), [1 3 2])).^2, 2));
    cnt = cnt + 1;
end

disp(' ')
disp('Time variants:')
arrayfun(@(x) disp([num2str(x) ': ' num2str(errors(x)), 'm/s']), 1:n_times);
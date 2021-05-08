%% Introduction
% This file computes the signal to noise ratio values for the chosen
% simulation environment.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Get parameters
params = generate_parameters();
load('./config/parameters/window.mat');
params.window = window;
params.output_folder = './images/snr';
time = compute_time(params);
load('./config/locations/test_sources.mat');
sources = test_sources;
clear test_sources window

%% Compute SNR
% Run in interations

% Determine iteration steps
memory = params.memory;
domain = params.domain;
sensors = params.sensors;
n_time_samples = domain.measurement_duration * sensors.sampling_rate;
one_sample = (3 * sensors.n_sensors * n_time_samples) * 8;
stepsize = floor(memory.sim_iteration_limit / one_sample);
n_sources = size(sources, 1);
starts = 1:stepsize:n_sources;
n_iterations = length(starts);

if params.print_messages
    disp('Start velocity generation and SNR computation')
    disp(['    using ', num2str(n_iterations), ' iterations for this setting']);
end

% Preallocate
snr.x = zeros(n_sources, 1, sensors.n_sensors);
snr.y = zeros(n_sources, 1, sensors.n_sensors);
snr.sources = zeros(n_sources, 3);

% Run the iterations
for idx = starts
    tic;
    indices = idx:min(idx+(stepsize-1), n_sources);
    
    velocity = compute_sensor_velocity(sources(indices,:), time, params);
    
    control = velocity;
    control.xt = zeros(size(control.xt));
    control.yt = zeros(size(control.yt));
    control = apply_sensor_noise(control, params);

    noisy_velocity = velocity;
    noisy_velocity.xt = velocity.xt + control.xt;
    noisy_velocity.yt = velocity.yt + control.yt;
    
    % Compute snr
    [snr.x(indices,:,:), snr.y(indices,:,:)] =...
        compute_snr(noisy_velocity, control, params);
    snr.sources(indices,:) = velocity.sources;
    if params.print_messages && idx == 1
        t = toc;
        disp(['    iteration in ', duration2str(t)]);
        disp(['        estimated time: ', duration2str(n_iterations * t)])
    end
end

%% Export to csv
cell_size = 0.02;
output_folder = './data/final_csv/';

snr.sources_bin = compute_bins(snr.sources, cell_size);

result = struct2table(snr);
result.Properties.VariableNames(1:2) = {'snr_vx' , 'snr_vy'};

result.x_source = result.sources(:, 1);
result.y_source = result.sources(:, 2);
result.orientation_source = result.sources(:, 3);
result.x_bin = result.sources_bin(:, 1);
result.y_bin = result.sources_bin(:, 2);
result.orientation_bin = result.sources_bin(:, 3);
result = removevars(result, {'sources', 'sources_bin'});

% Make it a tidy table
% - A column for the velocity component
result = stack(result, 1:2, 'NewDataVariableName', 'snr', 'IndexVariableName', 'velocity_component');
result.velocity_component = renamecats(result.velocity_component, {'x', 'y'});
% - A column for the sensor number
n_columns = width(result);
result = splitvars(result);
n_columns_after_split = width(result);
result = stack(result, n_columns:n_columns_after_split, 'NewDataVariableName', 'snr', 'IndexVariableName', 'sensor_number');
result.sensor_number = grp2idx(result.sensor_number);

writetable(result, [output_folder, 'snr.csv'])


%% Make overall plot
% conf.cell_size = 0.02;
% conf.color_range = [0, 60];
% conf.map = magma(diff(conf.color_range) / 5);
% conf.n_spokes = 7;
% conf.n_circles = 3;
% 
% snr.sources_bin = compute_bins(snr.sources, conf.cell_size);
% plot_snr(snr, params, conf);

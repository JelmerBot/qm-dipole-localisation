%% Introduction
% Generate training, validation and test source states using Poisson disc 
% sampling.

%% Setup
addpath(genpath('./functions'))
clear
clc

%% Get parameters
params = generate_parameters();
params.output_folder = './config/locations';
if ~exist(params.output_folder, 'dir')
    mkdir(params.output_folder)
end
save(fullfile(params.output_folder, 'params.mat'), 'params')

% High res sources
params.domain.min_source_distance = params.experiment.test_source_distance;
test_sources = compute_sampling_locations(params);
disp([num2str(params.experiment.test_source_distance) ' ' num2str(size(test_sources, 1))]);
save(fullfile(params.output_folder, 'test_sources.mat'), 'test_sources')

% Variable res sources
for idx = 1:length(params.experiment.train_source_distances)
    params.domain.min_source_distance = ...
        params.experiment.train_source_distances(idx);
    train_sources = compute_sampling_locations(params);
    disp([num2str(params.domain.min_source_distance) ' ' num2str(size(train_sources, 1))])
    save(fullfile(params.output_folder, ...
        ['train_sources_res', num2str(params.domain.min_source_distance),'.mat']), 'train_sources');
end

disp('done')
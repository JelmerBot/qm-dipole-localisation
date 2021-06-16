%% Introduction
% Evaluate all validated methods

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

params = generate_parameters();

%% Run data tests
for idx = 1:length(params.experiment.higher_noise_levels)
    test_methods_noise(idx)
end

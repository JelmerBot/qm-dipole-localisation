%% Introduction
% Evaluate all validated methods on the test sets. May take quite a
% while to complete!

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

params = generate_parameters();

%% Run data tests
for idx = 1:length(params.experiment.train_source_distances)
    test_methods_data(idx)
end

%% Run input tests
for idx = 1:length(params.experiment.input_modes)
    test_methods_input(idx)
end

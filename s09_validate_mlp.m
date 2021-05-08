%% Introduction
% Run the mlp crossvalidation for each trainingset. May take quite a
% while to complete for the larger training sets!

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

params = generate_parameters();

%% Validate KNN for training data
for idx = 1:length(params.experiment.train_source_distances)
    mlp_data_validation(idx);
end

%% Validate mlp for input modes
for idx = 1:length(params.experiment.input_modes)
    mlp_input_validation(idx);
end

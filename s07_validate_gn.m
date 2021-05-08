%% Introduction
% Run the GN crossvalidation for each trainingset. May take quite a
% while to complete for the larger training sets!

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Validate GN for training data

params = generate_parameters();
for idx = 1:length(params.experiment.train_source_distances)
    gn_data_validation(idx);
end

%% Validate GN for input modes

for idx = 1:length(params.experiment.input_modes)
    gn_input_validation(idx);
end

%% Introduction
% Run the NR crossvalidation for each trainingset. May take quite a
% while to complete for the larger training sets!

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Validate NR for training data

params = generate_parameters();
for idx = 1:length(params.experiment.train_source_distances)
    nr_data_validation(idx);
end

%% Validate NR for input modes

for idx = 1:length(params.experiment.input_modes)
    nr_input_validation(idx);
end

%% Introduction
% Computes and shows the average prediction and total training time

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

% Input folders (same values as in s11)
input_mode_folder = './data/20_05_21/';
data_folder = './data/20_05_21/';

% Load timing only of largest training set with x+y
load([input_mode_folder, 'timing_input_mode_1.mat'], 'timing');
input_mode_timing = struct2table(timing);
load([data_folder, 'timing_res0.01.mat'], 'timing');
data_timing = struct2table(timing);
clear timing;

% Load number of sources
load('config/locations/train_sources_res0.01.mat', 'train_sources')
load('config/locations/test_sources.mat', 'test_sources')
n_test_sources = size(test_sources, 1);
n_train_sources = size(train_sources, 1);
clear test_sources train_sources

% Time dimension was removed for all predictors except LCMV.
% Add that time to the total prediction time.
input_mode_timing{:,[6,7,10:end]} = input_mode_timing{:,[6,7,10:end]} +...
    input_mode_timing.reduction;
input_mode_timing.reduction = [];

% Add QM from the data_timing set
% Measurements could be from 2 different processors :(
input_mode_timing.qm = data_timing.reduction + data_timing.qm;

% Split training time and prediction time
prediction_times = input_mode_timing(:,6:end);
% average (methods were optimized for predicting large number of sources)
prediction_times{1,:} = prediction_times{1,:} ./ n_test_sources;

% Show avg prediction time and total train time in seconds
methods = prediction_times.Properties.VariableNames;
times = prediction_times{1,:};
prediction_times = table(methods', times',...
    'VariableNames', {'method', 'prediction_time'});
prediction_times = sortrows(prediction_times, 'prediction_time');
prediction_times.vertraging = prediction_times.prediction_time ./ prediction_times.prediction_time(2);
prediction_times 

train_times = input_mode_timing(:,1:4)

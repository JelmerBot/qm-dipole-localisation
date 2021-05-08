function mlp_data_validation(data_id)
%mlp_data_validation
%
%   Syntax mlp_data_validation(data_id)
%
% data_id refers to which trainingset to load.

narginchk(1, 1)
nargoutchk(0, 0)
if isa(data_id, 'char')
    data_id = round(str2double(data_id));
end

%% Parameters
params = generate_parameters();
% load('./config/parameters/window.mat', 'window');
params.window = @hamming; % Compiled code does not like the function loaded as above.

if data_id < 1 || data_id > length(params.experiment.train_source_distances)
    error('Invalid data_id');
end

sample_distance = params.experiment.train_source_distances(data_id);
input_file = ['./config/locations/train_sources_res',num2str(sample_distance) '.mat'];
load(input_file, 'train_sources');

output_file_results = ['./config/validation/mlp_res' num2str(sample_distance) '.mat'];
output_file_configuration = ['./config/parameters/mlp_res' num2str(sample_distance) '.mat'];
output_file_plot = ['./images/validation/mlp_res' num2str(sample_distance) '_'];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
% save('tmp_vel.mat', 'velocity');
% load('tmp_vel.mat', 'velocity');
velocity = apply_input_mode(velocity, params);
n_fold = params.experiment.n_fold;
n_samples = size(velocity.input, 1);
partitions = crossvalind('Kfold', n_samples, n_fold);

%% Search for best size
% random search
results = bayesopt_search(velocity, partitions, params);

results = extract_data(results);
save(output_file_results, 'results')

% load(output_file_results, 'results')

plot_variable_effects(results, output_file_plot, [false, true, true]);

[v, i] = min(results.mean_test_error);
mlp = params.mlp;
mlp.layer_sizes = repmat(results.layer_size(i), 1, results.n_layers(i));
mlp.learn_rate = results.learn_rate(i);
save(output_file_configuration, 'mlp');
disp(['Optimal learn rate: ', num2str(mlp.learn_rate)...
            ', num layers: ', num2str(results.n_layers(i)),...
            ', layer size: ', num2str(results.layer_size(i))])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
disp(' ')

end

function results = bayesopt_search(velocity, partitions, params)

mlp = params.mlp;

learn_rate = optimizableVariable('learn_rate', mlp.validate.learn_rate,...
    'Type', 'real', 'Transform', 'log');
n_layers = optimizableVariable('n_layers', mlp.validate.n_layers,...
    'Type', 'integer', 'Transform', 'none');
layer_size = optimizableVariable('layer_size', mlp.validate.layer_size,...
    'Type', 'integer', 'Transform', 'log');

fun = @(x) evaluate_mlp(x, @state, velocity, partitions, params);
results = bayesopt(fun,...
    [n_layers, layer_size, learn_rate],...
    'MaxObjectiveEvaluations', mlp.validate.evaluations,...
    'NumSeedPoints', ceil(mlp.validate.evaluations / 5),...
    'IsObjectiveDeterministic', true,...
    'PlotFcn', {});
end

function params = state(x, params)

params.mlp.learn_rate = x.learn_rate;
params.mlp.layer_sizes = repmat(x.layer_size, 1, x.n_layers);
params.mlp.verbose = false;

end

function [value, condition, user_data] = evaluate_mlp(x, fun, velocity, partitions, params)

% extract parameter values
params = fun(x, params);

% validate setting
info = cross_validate_mlp_setting(velocity, partitions, params);

% collect metrics
user_data.info = info;
user_data.mean_test_error = mean(info.test_errors(:, end));
user_data.std_test_error = std(info.test_errors(:,end));
user_data.mean_train_error = mean(info.train_errors(:, end));
user_data.std_train_error = std(info.train_errors(:, end));
total_train_times = sum(info.train_times, 2);
user_data.mean_train_time = mean(total_train_times);
user_data.std_train_time = std(total_train_times);

value = user_data.mean_test_error;
condition = [];
end

function results = extract_data(results)
% variables
names = results.XTrace.Properties.VariableNames;
txt = [];
for name = names
   eval([name{1} ' = results.XTrace.' name{1} ';'])
   txt = [txt name{1} ','];
end

% Metrics
udat = [results.UserDataTrace{:}];
mean_test_error = [udat.mean_test_error]';
std_test_error = [udat.std_test_error]';

mean_train_error = [udat.mean_train_error]';
std_train_error = [udat.std_train_error]';

mean_train_time = [udat.mean_train_time]';
std_train_time = [udat.std_train_time]';

% Output table
eval(['results = table(' txt ...
            'mean_test_error, std_test_error, '...
            'mean_train_error, std_train_error, '...
            'mean_train_time, std_train_time);']);
end

function plot_variable_effects(results, out_file, logtransform)
names = results.Properties.VariableNames;
n_variables = length(names) - 8;

for idx = 1:n_variables
    name = names{idx};
    name = strrep(name, '_', ' ');
    
    figure
    colors = get(gca, 'ColorOrder');
    subplot(1, 2, 1)
    plot(results{:,idx}, results.mean_test_error, '.', 'color', colors(2,:));
    hold on
    plot(results{:,idx}, results.mean_train_error, '.', 'color', colors(1,:));
    xlabel(name)
    ylabel('mean absolute error')
    if logtransform(idx)
        ax = gca;
        ax.XScale = 'log';
    end
    subplot(1, 2, 2)
    plot(results{:,idx}, results.mean_train_time, '.', 'color', colors(1,:))
    xlabel(name)
    ylabel('run time')
    if logtransform(idx)
        ax = gca;
        ax.XScale = 'log';
    end
    post_process_figure(1, 0.4, [0 1], [0.2 0.2])
    print([out_file, names{idx}, '.jpg'], '-djpeg', '-r600');
    close
end
end

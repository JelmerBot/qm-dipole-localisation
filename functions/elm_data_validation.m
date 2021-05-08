function elm_data_validation(data_id)
%elm_data_validation
%
%   Syntax elm_data_validation(data_id)
%
% data_id refers to which trainingset to load.

narginchk(1, 1)
nargoutchk(0, 0)
if isa(data_id, 'char')
    data_id = round(str2double(data_id));
end

%% Parameters
params = generate_parameters();
%load('./config/parameters/window.mat', 'window');
params.window = @hamming; % Compiled code does not like loading functions as above.

if data_id < 1 || data_id > length(params.experiment.train_source_distances)
    error('Invalid data_id');
end

sample_distance = params.experiment.train_source_distances(data_id);
input_file = ['./config/locations/train_sources_res',num2str(sample_distance) '.mat'];
load(input_file, 'train_sources');

output_file_results = ['./config/validation/elm_res' num2str(sample_distance) '.mat'];
output_file_configuration = ['./config/parameters/elm_res' num2str(sample_distance) '.mat'];
output_file_plot = ['./images/validation/elm_res' num2str(sample_distance) '.jpg'];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);
velocity = normalise_input(velocity);
[mean_absolute_error, run_time, n_nodes] = cross_validate_elm(velocity, params);
save(output_file_results, 'mean_absolute_error', 'run_time', 'n_nodes');

%% Analyse results

% load(output_file_results, 'mean_absolute_error', 'run_time', 'n_nodes');

figure
colors = get(gca, 'colororder');
subplot(1,2,1)
plot(n_nodes(1,:), mean(mean_absolute_error, 1), '-', 'color', colors(2,:));
hold on
plot(n_nodes, mean_absolute_error, '.', 'color', colors(1,:));
ax = gca; ax.XAxis.Scale = 'log';
xlabel('num nodes')
ylabel('mean absolute error')
ylim([0 2])
subplot(1,2,2)
plot(n_nodes(1,:), mean(run_time, 1), '-', 'color', colors(2,:));
hold on
plot(n_nodes, run_time, '.', 'color', colors(1,:));
ax = gca; ax.XAxis.Scale = 'log';
xlabel('num nodes')
ylabel('run time')
post_process_figure(1, 0.4, [0 0], [0 0])
print(output_file_plot, '-djpeg', '-r600');
close

upper_limit = find(any(isnan(mean_absolute_error), 1), 1, 'first');
if isempty(upper_limit)
    upper_limit = length(mean_absolute_error);
end
errors = mean(mean_absolute_error(1:upper_limit), 1);
[v, i] = min(errors);
elm = params.elm;
elm.n_nodes = n_nodes(1,i);
disp(['Optimal nodes: ', num2str(elm.n_nodes)])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
disp(' ')
save(output_file_configuration, 'elm');

end

function elm_input_validation(input_id)
%elm_input_validation
%
%   Syntax elm_input_validation(input_id)
%
% input_id refers to which input mode to use.

narginchk(1, 1)
nargoutchk(0, 0)
if isa(input_id, 'char')
    input_id = round(str2double(input_id));
end

%% Parameters
params = generate_parameters();
% load('./config/parameters/window.mat', 'window');
params.window = @hamming; % Compiled code does not like loading functions as above.
if input_id < 1 || input_id > length(params.experiment.input_modes)
    error('Invalid input_id');
end
params.sensors.input_mode = params.experiment.input_modes{input_id};

sample_distance = min(params.experiment.train_source_distances);
input_file = ['./config/locations/train_sources_res',num2str(sample_distance) '.mat'];
load(input_file, 'train_sources');

output_file_results = ['./config/validation/elm_input_mode_' num2str(input_id) '.mat'];
output_file_configuration = ['./config/parameters/elm_input_mode_' num2str(input_id) '.mat'];
output_file_plot = ['./images/validation/elm_input_mode_' num2str(input_id)];

% %% Velocity and validation
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

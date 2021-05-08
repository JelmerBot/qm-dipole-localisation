function knn_data_validation(data_id)
%knn_data_validation
%
%   Syntax knn_data_validation(data_id)
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
params.window = @hamming; % Compiled code does not like loading functions as above.

if data_id < 1 || data_id > length(params.experiment.train_source_distances)
    error('Invalid data_id');
end

sample_distance = params.experiment.train_source_distances(data_id);
input_file = ['./config/locations/train_sources_res',num2str(sample_distance) '.mat'];
load(input_file, 'train_sources');

output_file_results = ['./config/validation/knn_res' num2str(sample_distance) '.mat'];
output_file_configuration = ['./config/parameters/knn_res' num2str(sample_distance) '.mat'];
output_file_plot = ['./images/validation/knn_res' num2str(sample_distance) '.jpg'];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);
velocity = normalise_input(velocity);
[mean_absolute_error, run_time, k] = cross_validate_knn(velocity, params);
save(output_file_results, 'mean_absolute_error', 'run_time', 'k');

%% Analyse results

% load(output_file_results, 'mean_absolute_error', 'run_time', 'k');

figure
colors = get(gca, 'colororder');
subplot(1,2,1)
plot(k(1,:), mean(mean_absolute_error, 1), '-', 'color', colors(2,:));
hold on
plot(k, mean_absolute_error, '.', 'color', colors(1,:));
xlabel('k')
ylabel('mean absolute error')
ylim([0 2])
subplot(1,2,2)
plot(k(1,:), mean(run_time, 1), '-', 'color', colors(2,:));
hold on
plot(k, run_time, '.', 'color', colors(1,:));
xlabel('k')
ylabel('run time')
post_process_figure(1, 0.4, [0 0], [0 0])
print(output_file_plot, '-djpeg', '-r600');
close

errors = mean(mean_absolute_error, 1);
[v, i] = min(errors);
knn = params.knn;
knn.k = k(1,i);
disp(['Optimal k: ', num2str(knn.k)])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
disp(' ')
save(output_file_configuration, 'knn');

end

function cwt_data_validation(data_id)
%cwt_data_validation
%
%   Syntax cwt_data_validation(data_id)
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
disp(['Num training_sources: ', num2str(size(train_sources, 1))])

output_file_results = ['./config/validation/cwt_res' num2str(sample_distance) '.mat'];
output_file_configuration = ['./config/parameters/cwt_res' num2str(sample_distance) '.mat'];
output_file_plot = ['./images/validation/cwt_res' num2str(sample_distance)];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
% save('tmp_vel.mat', 'velocity')
% load('tmp_vel.mat', 'velocity')
wavelets = compute_wavelets(train_sources, params);

[mean_absolute_error, run_time, t_min, t_max, c_x, c_y] = cross_validate_cwt(velocity, wavelets, params);
save(output_file_results, 'mean_absolute_error', 'run_time', 't_min', 't_max', 'c_x', 'c_y');

%% Analyse results

pre_process_figure

% load(output_file_results, 'mean_absolute_error', 'run_time', 't_min', 't_max', 'c_x', 'c_y');

figure
colors = get(gca, 'colororder');
subplot(1,2,1)
plot3(t_min, t_max, mean_absolute_error, '.', 'color', colors(1, :), 'markersize', 3);
hold on
plot3(t_min(:, 1), t_max(:, 1), mean(mean_absolute_error, 2), '.', 'color', colors(2, :));
grid on
view(-55.9, 10)
xlabel('$t_{min}$')
ylabel('$t_{max}$')
zlabel('mean absolute error')
subplot(1,2,2)
plot3(t_min, t_max, run_time, '.', 'color', colors(1, :), 'markersize', 3);
hold on
plot3(t_min(:, 1), t_max(:, 1), mean(run_time, 2), '.', 'color', colors(2, :));
grid on
view(-55.9, 10)
xlabel('$t_{min}$')
ylabel('$t_{max}$')
zlabel('run time')
post_process_figure(1, 0.4, [0 0], [0 0])
print([output_file_plot '.jpg'], '-djpeg', '-r600');
close

names = {'t_min', 't_max', 'c_x', 'c_y'};
results = [t_min(:, 1), t_max(:, 1), c_x(:, 1), c_y(:, 1)];
for idx = 1:length(names)
    name = names{idx};
    name = strrep(name, '_', ' ');
    
    figure
    colors = get(gca, 'ColorOrder');
    subplot(1, 2, 1)
    plot(results(:, idx), mean_absolute_error, '.', 'color', colors(1,:), 'markersize', 3);
    hold on
    plot(results(:, idx), mean(mean_absolute_error, 2), '.', 'color', colors(2,:));
    xlabel(name)
    ylabel('mean absolute error')

    subplot(1, 2, 2)
    plot(results(:, idx), run_time, '.', 'color', colors(1,:), 'markersize', 3);
    hold on
    plot(results(:, idx), mean(run_time, 2), '.', 'color', colors(2,:))
    xlabel(name)
    ylabel('run time')

    post_process_figure(1, 0.4, [0 1], [0.2 0.2])
    print([output_file_plot, '_', names{idx}, '.jpg'], '-djpeg', '-r600');
    close
end

errors = mean(mean_absolute_error, 2);
[v, i] = min(errors);
cwt = params.cwt;
cwt.threshold_min = t_min(i, 1);
cwt.threshold_max = t_max(i, 1);
cwt.c_x = c_x(i, 1);
cwt.c_y = c_y(i, 1);
disp(['Optimal tmin: ', num2str(cwt.threshold_min)...
            ', tmax: ', num2str(cwt.threshold_max)...
            ', cx: ', num2str(cwt.c_x)...
            ', cy: ', num2str(cwt.c_y)])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
disp(' ')
save(output_file_configuration, 'cwt');

end

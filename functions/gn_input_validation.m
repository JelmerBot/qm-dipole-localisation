function gn_input_validation(input_id)
%gn_input_validation
%
%   Syntax gn_input_validation(input_id)
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

output_file_results = ['./config/validation/gn_input_mode_' num2str(input_id) '.mat'];
output_file_configuration = ['./config/parameters/gn_input_mode_' num2str(input_id) '.mat'];
output_file_plot = ['./images/validation/gn_input_mode_' num2str(input_id) '.jpg'];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);
[mean_absolute_error, run_time, y] = validate_gn(velocity, params);
save(output_file_results, 'mean_absolute_error', 'run_time', 'y');

%% Analyse results

% load(output_file_results, 'mean_absolute_error', 'run_time', 'y');

figure
subplot(1,2,1)
plot(y(:), mean_absolute_error(:), '.');
xlabel('y')
ylabel('mean absolute error')
ylim([0 2])
subplot(1,2,2)
plot(y(:), run_time(:), '.');
xlabel('y')
ylabel('run time')
post_process_figure(1, 0.4, [0 0], [0 0])
print(output_file_plot, '-djpeg', '-r600');
close

[v, i] = min(mean_absolute_error);
gn = params.gn;
gn.starting_estimate(2) = y(i);
disp(['starting y: ', num2str(gn.starting_estimate(2))])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
save(output_file_configuration, 'gn');
disp(' ')

end

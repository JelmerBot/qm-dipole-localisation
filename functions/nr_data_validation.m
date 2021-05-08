function nr_data_validation(data_id)
%nr_data_validation
%
%   Syntax nr_data_validation(data_id)
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

output_file_results = ['./config/validation/nr_res' num2str(sample_distance) '.mat'];
output_file_configuration = ['./config/parameters/nr_res' num2str(sample_distance) '.mat'];
output_file_plot = ['./images/validation/nr_res' num2str(sample_distance) '_'];

%% Velocity and validation
velocity = generate_noisy_reduced_velocity(train_sources, params);
velocity = apply_input_mode(velocity, params);

results = validate_nr(velocity, params);

mean_absolute_error = results.ObjectiveTrace;
starting_y = results.XTrace.starting_y;
% step_size = results.XTrace.step_size;
norm_limit = results.XTrace.norm_limit;
% regularization = results.XTrace.regularization;

udat = [results.UserDataTrace{:}];
run_time = [udat.run_time]';

results = table(starting_y, norm_limit, mean_absolute_error, run_time);

save(output_file_results, 'results');

%% Save configuration

% load(output_file_results, 'results');

[v, i] = min(results.mean_absolute_error);
nr = params.nr;
nr.starting_estimate(2) = results.starting_y(i);
nr.norm_limit = results.norm_limit(i);
save(output_file_configuration, 'nr');
disp(' ')
disp(['Optimal y: ', num2str(nr.starting_estimate(2)) ', l: ', num2str(nr.norm_limit)])
disp([9, 'mean absolute error:' 9 9 num2str(v)])
disp(' ')

%% Analyse results

logtransform = [false false];
names = results.Properties.VariableNames;
n_variables = length(names) - 2;

for idx = 1:n_variables
    name = names{idx};
    name = strrep(name, '_', ' ');
    
    figure
    subplot(1, 2, 1)
    plot(results{:,idx}, results.mean_absolute_error, '.')
    xlabel(name)
    ylabel('mean absolute error')
    if logtransform(idx)
        ax = gca;
        ax.XScale = 'log';
    end
    subplot(1, 2, 2)
    plot(results{:,idx}, results.run_time, '.')
    xlabel(name)
    ylabel('run time')
    if logtransform(idx)
        ax = gca;
        ax.XScale = 'log';
    end
    post_process_figure(1, 0.4, [0 0], [0 0])
    print([output_file_plot, names{idx} '.jpg'], '-djpeg', '-r600');
    close
end

end

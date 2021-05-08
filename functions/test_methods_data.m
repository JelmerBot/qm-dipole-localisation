function test_methods_data(data_id)
%test_methods_data
%
%   Syntax test_methods_data(data_id)
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

if ~exist(params.output_folder, 'dir')
    mkdir(params.output_folder)
end

% Load sources
sample_distance = params.experiment.train_source_distances(data_id);
input_file = ['./config/locations/train_sources_res',num2str(sample_distance) '.mat'];
load(input_file, 'train_sources');
input_file = './config/locations/test_sources.mat';
load(input_file, 'test_sources');

% Load validated parameters
loading_name = @(x) ['./config/parameters/' x '_res' num2str(sample_distance) '.mat'];
load(loading_name('cwt'), 'cwt');
params.cwt = cwt;
load(loading_name('elm'), 'elm');
params.elm = elm;
load(loading_name('gn'), 'gn');
params.gn = gn;
load(loading_name('knn'), 'knn');
params.knn = knn;
load(loading_name('nr'), 'nr');
params.nr = nr;
load(loading_name('mlp'), 'mlp');
params.mlp = mlp;
params.mlp.verbose = false;

save(fullfile(params.output_folder,...
              ['params_res' num2str(sample_distance) '.mat']),...
     'params');

clear cwt elm gn knn nr mlp

%% Generate training data and train ELM | MLP

show('Generate training velocity', params)
tic
train_velocity = generate_noisy_reduced_velocity(train_sources, params);
train_velocity = apply_input_mode(train_velocity, params);
timing.velocity_train = toc;
tic
wavelets = compute_wavelets(train_sources, params);
timing.wavelets = toc;
show('    done', params)

show('Training elm', params)
tic
elm = train_elm(normalise_input(train_velocity), params);
timing.elm_train = toc;
show('    done', params)

% Training the mlp requires validation set (80-20 split).
show('Training mlp', params)
partitions = crossvalind('Kfold', size(train_velocity.sources, 1), params.experiment.n_fold);
validation_idx = partitions == 2;
train_idx = ~validation_idx;
mlp_validation_velocity = select_velocity_indices(train_velocity, validation_idx);
mlp_train_velocity = select_velocity_indices(train_velocity, train_idx);
tic
mlp = train_mlp(mlp_validation_velocity, mlp_train_velocity, params);
timing.mlp_train = toc;
show('    done', params)

%% Loop over testing locations

% determine loop size
memory = params.memory;
domain = params.domain;
sensors = params.sensors;
n_time_samples = domain.measurement_duration * sensors.sampling_rate;
one_sample = (3 * sensors.n_sensors * n_time_samples) * 8;
stepsize = floor(memory.sim_iteration_limit / one_sample);
n_sources = size(test_sources, 1);
starts = 1:stepsize:n_sources;
n_iterations = length(starts);

if params.print_messages
    disp('Start prediction loop')
    disp(['    using ', num2str(n_iterations), ' iterations']);
end

time = compute_time(params);
gcp;

cnt = 0;
for idx = starts
    cnt = cnt + 1;
    t1 = tic();
    show(['Iteration ' num2str(cnt)], params);
    
    indices = idx:min(idx+(stepsize-1), n_sources);
    n_indices = length(indices);
    
    % Compute test velocity
    show('    computing velocity', params)
    tic
    velocity = compute_sensor_velocity(test_sources(indices,:), time, params);
    velocity = apply_sensor_noise(velocity, params);
    timing.velocity_test = toc;
    
    % Compute lcmv predictions
    show('    computing lcmv predictions', params)
    tic
    train = normalise_input(apply_input_mode(train_velocity, params));
    test = apply_input_mode(velocity, params);
    lcmv_predictions = predict_lcmv(test, train, params);
    timing.lcmv = toc;
    
    % Extract signal (removes time dimension)
    show('    removing time dimension', params)
    tic;
    velocity = extract_signal_over_sensors(velocity, params);
    velocity = apply_input_mode(velocity, params);
    timing.reduction = toc;

    % Compute cwt predictions
    show('    computing cwt predictions', params)
    tic
    cwt_predictions = predict_cwt(velocity, wavelets, params);
    timing.cwt = toc;
    
    % Compute gn predictions
    show('    computing gn predictions', params)
    tic
    gn_predictions = predict_gn(velocity, params);
    timing.gn = toc;
    
    % Compute nr predictions
    show('    computing nr predictions', params)
    tic
    nr_predictions = predict_nr(velocity, params);
    timing.nr = toc;
    
    % Compute quadrature predictions
    show('    computing qm predictions', params)
    tic
    quad_predictions = predict_qm(velocity, params);
    timing.qm = toc;
    
    % Compute random predictions
    show('    computing random predictions', params)
    tic
    rand_predictions = predict_random(velocity, params);
    timing.random = toc;
    
    % Compute elm predictions
    show('    computing elm predictions', params)
    tic;
    elm_predictions = predict_elm(normalise_input(velocity), elm, params);
    timing.elm_predict = toc;
    
    % Compute mlp predictions
    show('    computing mlp predictions', params)
    tic;
    mlp_predictions = predict_mlp(velocity, mlp, params);
    timing.mlp_predict = toc;

    % Compute knn predictions
    show('    computing knn predictions', params)
    tic;
    knn_predictions = predict_knn(normalise_input(velocity), normalise_input(train_velocity), params);
    timing.knn = toc;
       
    % Write predictions to file
    show('    saving predictions to file', params)
    method = {'LCMV', 'CWT', 'GN', 'NR', 'QM', 'RND', 'ELM', 'MLP', 'KNN'};
    n_methods = length(method);
    method = repelem(method, 1, n_indices);
    source = repmat(velocity.sources, n_methods, 1);
    prediction = [lcmv_predictions; cwt_predictions; gn_predictions;...
                  nr_predictions; quad_predictions; rand_predictions;...
                  elm_predictions; mlp_predictions; knn_predictions];
    predictions = table(method', source, prediction);
    save(fullfile(params.output_folder,...
                  ['predictions_res' num2str(sample_distance) '_' num2str(cnt) '.mat']),...
         'predictions', '-v7.3');
    save(fullfile(params.output_folder,...
                  ['timing_res' num2str(sample_distance) '_' num2str(cnt) '.mat']),...
         'timing', '-v7.3');
                    
    if params.print_messages && cnt == 1
        t = toc(t1);
        disp(['    iteration in ', duration2str(t)]);
        disp(['        estimated time: ', duration2str(n_iterations * t)])
    end
end

show('Merging iteration files...', params)
% Merge iteration files
merged = [];
times = [];
for idx = 1:n_iterations
    name = fullfile(params.output_folder,...
                   ['predictions_res' num2str(sample_distance) '_' num2str(idx) '.mat']);
    load(name, 'predictions')
    delete(name);
    name = fullfile(params.output_folder,...
                   ['timing_res' num2str(sample_distance) '_' num2str(idx) '.mat']);
    load(name, 'timing')
    delete(name);

    merged = [merged; predictions];
    times = [times, timing];
end
clear timing;
timing.velocity_train = times(1).velocity_train;
timing.wavelets = times(1).wavelets;
timing.elm_train = times(1).elm_train;
timing.mlp_train = times(1).mlp_train;
timing.velocity_test = sum([times.velocity_test]);
timing.elm_predict = sum([times.elm_predict]);
timing.mlp_predict = sum([times.mlp_predict]);
timing.lcmv = sum([times.lcmv]);
timing.reduction = sum([times.reduction]);
timing.cwt = sum([times.cwt]);
timing.gn = sum([times.gn]);
timing.nr = sum([times.nr]);
timing.qm = sum([times.qm]);
timing.knn = sum([times.knn]);
timing.random = sum([times.random]);

save(fullfile(params.output_folder,...
              ['timing_res' num2str(sample_distance) '.mat']),...
     'timing');

predictions = merged;
clear merged;
save(fullfile(params.output_folder,...
              ['predictions_res' num2str(sample_distance) '.mat']),...
     'predictions');

 
show('Done!', params)
 
end

function show(msg, params)
    if params.print_messages
       disp(msg); 
    end
end

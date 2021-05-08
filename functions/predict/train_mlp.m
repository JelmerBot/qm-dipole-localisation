function [mlp, info] = train_mlp(validate_velocity, training_velocity, params)
%train_mlp
%
% Syntax: [mlp, info] = train_mlp(validate_velocity, training_velocity, params)
%
% Implements MLP training.
narginchk(3, 3)
nargoutchk(1, 2)
validate_velocity.time = 0; % force that time is removed
training_velocity.time = 0; % force that time is removed
validatevelocity(validate_velocity)
validatevelocity(training_velocity)

% Extract parameter values
mlp = params.mlp;
levels = mlp.training_levels;
n_levels = length(levels);

% Apply velocity normalisation
if mlp.input_normalisation
   validate_velocity = normalise_input(validate_velocity);
   training_velocity = normalise_input(training_velocity);
end

% Configure MLP
n_inputs = size(validate_velocity.input, 3);
n_outputs = size(validate_velocity.sources, 2);
layers = configure_mlp(n_inputs, n_outputs, params);

% Preallocate values
info.test_errors = zeros(1, n_levels);
info.train_errors = zeros(1, n_levels);
info.train_times = zeros(1, n_levels);
info.training_info = cell(1, n_levels);
info.networks = cell(1, n_levels);

% Run staged training
for idx = 1:n_levels
    lvl = levels(idx); 
    
    if mlp.verbose
        disp(['Training level ' num2str(idx) ': ' num2str(lvl)])
    end
    
    % Filter out sources further than lvl from the origin
    training_mask = sqrt(sum(training_velocity.sources(:, 1:2).^2, 2)) < lvl;
    train_level_velocity = select_velocity_indices(training_velocity, training_mask);
    validate_mask = sqrt(sum(validate_velocity.sources(:, 1:2).^2, 2)) < lvl;
    validate_level_velocity = select_velocity_indices(validate_velocity, validate_mask);
    
    % Reformat for training
    train_input = permute(train_level_velocity.input, [3 2 4 1]);
    train_output = train_level_velocity.sources;
    validate_input = permute(validate_level_velocity.input(:,:,:), [3 2 4 1]);
    validate_output = validate_level_velocity.sources;
    
    % Train the network
    t = tic;
    options = configure_training({validate_input, validate_output}, params);
    [trained_net, info.training_info{idx}] = trainNetwork(train_input, train_output, layers, options);
    info.train_times(idx) = toc(t);
    layers = trained_net.Layers;
    info.networks{idx} = trained_net;
    
    % Compute training and testing errors
    predictions = predict_mlp(validate_velocity, trained_net, params);
    info.test_errors(idx) = mean(compute_prediction_errors(predictions, validate_velocity.sources, params));
    
    predictions = predict_mlp(training_velocity, trained_net, params);
    info.train_errors(idx) = mean(compute_prediction_errors(predictions, training_velocity.sources, params));
end

mlp = trained_net;

end

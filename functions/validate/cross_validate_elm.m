function [mean_absolute_error, run_time, n_nodes] = cross_validate_elm(velocity, params)
%cross_validate_knn
%
% Syntax: [mean_absolute_error, run_time, n_nodes] = cross_validate_elm(velocity, params)
narginchk(2, 2)
nargoutchk(3, 3)
validatevelocity(velocity)

n_values = params.elm.validate.evaluations;
n_fold = params.experiment.n_fold;
n_samples = size(velocity.input, 1);

% Preallocate data
mean_absolute_error = zeros(n_fold, n_values);
run_time = zeros(n_fold, n_values);
n_nodes = zeros(n_fold, n_values);

% Partition the data
partitions = crossvalind('Kfold', n_samples, n_fold);

for idx = 1:n_fold
    validate_idx = partitions == idx;
    train_idx = ~validate_idx;
    
    train = select_velocity_indices(velocity, train_idx);
    validate = select_velocity_indices(velocity, validate_idx);
    
    [...
        mean_absolute_error(idx, :),...
        run_time(idx, :),...
        n_nodes(idx, :)...
    ] = validate_elm(validate, train, params);
end

end
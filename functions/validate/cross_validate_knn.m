function [mean_absolute_error, run_time, k] = cross_validate_knn(velocity, params)
%cross_validate_knn
%
% Syntax: [mean_absolute_error, run_time, k] = cross_validate_knn(velocity, params)
narginchk(2, 2)
nargoutchk(3, 3)
validatevelocity(velocity)

n_fold = params.experiment.n_fold;
n_k = length(params.knn.validate.k);
n_samples = size(velocity.input, 1);

% Preallocate data
mean_absolute_error = zeros(n_fold, n_k);
run_time = zeros(n_fold, n_k);
k = params.knn.validate.k;

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
        k(idx, :)...
    ] = validate_knn(validate, train, params);
end

end
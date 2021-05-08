function [training_info] = cross_validate_mlp_setting(velocity, partitions, params)
%cross_validate_mlp_setting
%
% Syntax: [training_info] = cross_validate_mlp_setting(velocity, partitions, params)
%
% Trains and validates an mlp using the provided parameters
% 5-fold 80-20 crossvalidation is used
narginchk(3, 3)
nargoutchk(1, 1)
validatevelocity(velocity)

% Partition the data
n_fold = params.experiment.n_fold;
n_levels = length(params.mlp.training_levels);

% Preallocate results
training_info.test_errors = zeros(n_fold, n_levels);
training_info.train_errors = zeros(n_fold, n_levels);
training_info.train_times = zeros(n_fold, n_levels);
training_info.training_info = cell(n_fold, n_levels); % For debugging

for idx = 1:n_fold
    
    if params.mlp.verbose
        disp(' ')
        disp(['Validation round: ' num2str(idx)])
    end
    
    validate_idx = partitions == idx;
    train_idx = ~validate_idx;
    
    train = select_velocity_indices(velocity, train_idx);
    validate = select_velocity_indices(velocity, validate_idx);

    [~, info] = train_mlp(validate, train, params);
    
    training_info.train_times(idx, :) = info.train_times;
    training_info.test_errors(idx, :) = info.test_errors;
    training_info.train_errors(idx, :) = info.train_errors;
    training_info.training_info(idx, :) = info.training_info;
end



end
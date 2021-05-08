function [mean_absolute_error, run_time, t_min, t_max, c_x, c_y] = cross_validate_cwt(velocity, wavelets, params)
%cross_validate_cwt
%
% Syntax: [mean_absolute_error, run_time, t_min, t_max, c_x, c_y] = cross_validate_cwt(velocity, wavelets, params)
narginchk(3, 3)
nargoutchk(6, 6)
validatevelocity(velocity)


% Perform parameter search
gcp;
cwt = params.cwt;

t_min = optimizableVariable('t_min', cwt.validate.threshold_min,...
            'Type', 'real', 'Transform', 'none');
t_max = optimizableVariable('t_max', cwt.validate.threshold_max,...
            'Type', 'real', 'Transform', 'none');
c_x = optimizableVariable('c_x', cwt.validate.c_x,...
            'Type', 'real', 'Transform', 'none');
c_y = optimizableVariable('c_y', cwt.validate.c_y,...
            'Type', 'real', 'Transform', 'none');   
        
fun = @(x) evaluate_cwt(x, velocity, wavelets, params);
results = bayesopt(fun,...
    [t_min, t_max, c_x, c_y],...
    'MaxObjectiveEvaluations', params.cwt.validate.evaluations,...
    'IsObjectiveDeterministic', true,...
    'NumCoupledConstraints', 1,...
    'PlotFcn', {});

% Extract usefull data
udat = [results.UserDataTrace{:}];
mean_absolute_error = [udat.mean_absolute_errors]';
run_time = [udat.run_times]';
t_min = repmat(results.XTrace.t_min, 1, 5);
t_max = repmat(results.XTrace.t_max, 1, 5);
c_x = repmat(results.XTrace.c_x, 1, 5);
c_y = repmat(results.XTrace.c_y, 1, 5);

end

function [value, condition, user_data] = evaluate_cwt(x, velocity, wavelets, params)

n_samples = size(velocity.xt, 1);
n_fold = params.experiment.n_fold;

% t_min has to be lower than t_max
condition = x.t_min - x.t_max;
if condition >= 0
    value = inf;
    user_data.run_times = inf(n_fold, 1);
    user_data.mean_absolute_errors = inf(n_fold, 1);
    return
end

% set values
params.cwt.threshold_min = x.t_min;
params.cwt.threshold_max = x.t_max;
params.cwt.c_x = x.c_x;
params.cwt.c_y = x.c_y;

partitions = crossvalind('Kfold', n_samples, n_fold);
mean_absolute_errors = zeros(n_fold, 1);
run_times = zeros(n_fold, 1);

% Perform cross vallidation
for idx = 1:n_fold
    validate_idx = partitions == idx;
    train_idx = ~validate_idx;
    
    wavs = select_wavelet_indices(wavelets, train_idx);
    validate = select_velocity_indices(velocity, validate_idx);
    
    t = tic();
    predictions = predict_cwt(validate, wavs, params);
    run_times(idx) = toc(t);
    
    mean_absolute_errors(idx) = mean(compute_prediction_errors(predictions, validate.sources, params));    
end

% Put data in BayesOpt format
user_data.run_times = run_times;
user_data.mean_absolute_errors = mean_absolute_errors;

value = mean(mean_absolute_errors);
end
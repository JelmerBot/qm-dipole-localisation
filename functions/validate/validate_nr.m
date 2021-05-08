function results  = validate_nr(validate, params)
%validate_nr
%
% Syntax: results  = validate_nr(validate, params)
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(validate)

gcp;
nr = params.nr;

starting_y = optimizableVariable('starting_y', nr.validate.starting_y,...
    'Type', 'real', 'Transform', 'none');
% step_size = optimizableVariable('step_size', nr.validate.step_size,...
%     'Type', 'real', 'Transform', 'log');
norm_limit = optimizableVariable('norm_limit', nr.validate.norm_limit,...
    'Type', 'real', 'Transform', 'log');
% regularization = optimizableVariable('regularization', nr.validate.regularization,...
%     'Type', 'real', 'Transform', 'log');

fun = @(x) evaluate_nr(x, validate, params);
results = bayesopt(fun,...
    [starting_y, norm_limit],...
    'MaxObjectiveEvaluations', nr.validate.evaluations,...
    'NumSeedPoints', ceil(nr.validate.evaluations / 5),... 
    'IsObjectiveDeterministic', true,...
    'PlotFcn', {});
end

function [value, condition, user_data] = evaluate_nr(x, validate, params)

% extract parameter values
params.nr.starting_estimate(2) = x.starting_y;
% params.nr.step_size = x.step_size;
params.nr.norm_limit = x.norm_limit;
% params.nr.regularization = x.regularization;

% compute predictions
t = tic();
predictions = predict_nr(validate, params);
user_data.run_time = toc(t);

% Evaluate predictions
value = mean(compute_prediction_errors(predictions, validate.sources, params));
condition = [];
end

function [mean_absolute_error, run_time, y]  = validate_gn(validate, params)
%validate_gn
%
% Syntax: [mean_absolute_error, run_time, y]  = validate_gn(validate, params)
narginchk(2, 2)
nargoutchk(3, 3)
validatevelocity(validate)

gcp;
gn = params.gn;

% Determine y-values
y = linspace(min(gn.validate.starting_y), max(gn.validate.starting_y), gn.validate.evaluations);
n_y = length(y);
mean_absolute_error = zeros(1, n_y);
run_time = zeros(1, n_y);

disp(' ')
disp('Validate GN')
disp(['y in [', num2str(min(y))...
            ' ' num2str(max(y)) ']']);
disp(' ')
disp(['y', 9, 'mean absolute error', 9, 'run time'])
for idx = 1:n_y
    y_iter = y(idx);
    [...
        mean_absolute_error(idx),...
        run_time(idx)...
    ] = evaluate_gn(y_iter, validate, params);
    disp([num2str(y_iter), 9,...
          num2str(round(mean_absolute_error(idx), 4)), 9,9,...
          duration2str(round(run_time(idx), 4))])
end
disp(' ')
disp('done')
disp(' ')

end

function [mean_absolute_error, run_time] = evaluate_gn(y, validate, params)

params.gn.starting_estimate(2) = y;

t = tic();
predictions = predict_gn(validate, params);
run_time = toc(t);

mean_absolute_error = mean(compute_prediction_errors(predictions, validate.sources, params));

end
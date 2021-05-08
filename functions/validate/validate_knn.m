function [mean_absolute_error, run_time, k]  = validate_knn(validate, train, params)
%validate_knn
%
% Syntax: [mean_absolute_error, run_time, k]  = validate_knn(validate, train, params)
narginchk(3, 3)
nargoutchk(3, 3)
validatevelocity(validate)
validatevelocity(train)

gcp;
n_k = length(params.knn.validate.k < size(validate.input, 1));
mean_absolute_error = zeros(1, n_k);
run_time = zeros(1, n_k);
k = params.knn.validate.k;

disp(' ')
disp('Validate KNN')
disp(['k in [', num2str(min(params.knn.validate.k))...
            ' ' num2str(max(params.knn.validate.k)) ']']);
disp(' ')
disp(['k', 9, 'mean absolute error', 9, 'run time'])
for idx = 1:n_k
    k_iter = k(idx);
    [...
        mean_absolute_error(idx),...
        run_time(idx)...
    ] = evaluate_k(k_iter, validate, train, params);
    disp([num2str(k_iter), 9,...
          num2str(round(mean_absolute_error(idx), 4)), 9,9,...
          duration2str(round(run_time(idx), 4))])
end
disp(' ')
disp('done')
disp(' ')

end

function [mean_absolute_error, run_time] = evaluate_k(k, validate, train, params)

params.knn.k = k;
t = tic();
predictions = predict_knn(validate, train, params);
run_time = toc(t);

mean_absolute_error = mean(compute_prediction_errors(predictions, validate.sources, params));

end
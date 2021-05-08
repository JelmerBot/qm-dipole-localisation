function [mean_absolute_error, run_time, n_nodes] = validate_elm(validate, train, params)
%validate_knn
%
% Syntax: [mean_absolute_error, run_time, n_nodes] = validate_knn(validate, train, params)
narginchk(3, 3)
nargoutchk(3, 3)
validatevelocity(validate)
validatevelocity(train)

% Compute logarithmic spaces n_nodes
n_sources = size(train.input, 1) + size(validate.input, 1);
max_nodes_in_memory = floor(sqrt(params.memory.elm_nodes_limit / 8));
nodes_limit = min(n_sources, max_nodes_in_memory);
exponents = linspace(1, log10(nodes_limit), params.elm.validate.evaluations);
n_nodes = round(10.^exponents); % some duplicates may exist not really an issue though

% preallocate results
mean_absolute_error = nan(1, length(n_nodes));
run_time = nan(1, length(n_nodes));

% Reset warnings
lastwarn('');

disp(' ')
disp('Validate ELM')
disp(['n in [', num2str(min(n_nodes)) ' ' num2str(max(n_nodes)) ']']);
disp(' ')
disp(['n', 9, 'mean absolute error', 9, 'run time'])

% Evaluation loop
for idx = 1:length(n_nodes)
    nodes = n_nodes(idx);
    [mean_absolute_error(idx), run_time(idx)] = evaluate_elm(nodes, validate, train, params);
    disp([num2str(nodes), 9,...
          num2str(round(mean_absolute_error(idx), 4)), 9,9,...
          duration2str(round(run_time(idx), 4))])
    % Early stopping on error
    if mean_absolute_error(idx) > 2
        break;
    end
    [~, id] = lastwarn;
    if strcmp(id, 'MATLAB:nearlySingularMatrix') % || strcmp(id, 'ELM:notEnoughSources')
        break;
    end
end
disp(' ')
disp('done')
disp(' ')

end

function [mean_absolute_error, run_time] = evaluate_elm(n_nodes, validate, train, params)

params.elm.n_nodes = n_nodes;
t = tic();
predictions = train_and_predict_elm(validate, train, params);
run_time = toc(t);

mean_absolute_error = mean(compute_prediction_errors(predictions, validate.sources, params));

end
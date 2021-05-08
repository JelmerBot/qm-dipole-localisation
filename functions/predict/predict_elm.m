function predictions = predict_elm(testing_velocity, elm, params)
%predict_elm
%
% Syntax: predictions = predict_elm(testing_velocity, elm, params)
%
% Implements ELM predictions.
narginchk(3, 3)
nargoutchk(1, 1)
testing_velocity.time = 0;
validatevelocity(testing_velocity)

data = permute(testing_velocity.input, [1 3 2]);
validateattributes(data, {'double'}, {'ncols', elm.n_inputs}, 'predict', 'data', 1);
n_sources = size(data, 1);
          
% Preallocate output
predictions = zeros(n_sources, elm.n_outputs);

% Iterations to prevent memory issues
one_sample_size = elm.n_nodes * 8;
step_size = floor(params.memory.elm_iteration_limit / one_sample_size);
starts = 1:step_size:n_sources;

% Prediction loop
for start = starts
    indices = start:min(start+(step_size-1), n_sources);

    V = data(indices, :);
    H = elm.f_internal((elm.weights_internal * V')' + elm.bias_internal);
    predictions(indices, :) = H*elm.weights_output;
end

end
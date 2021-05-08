function elm = train_elm(training_velocity, params)
%train_elm
%
% Syntax: elm = train_elm(training_velocity, params)
%
% Implements ELM training.
narginchk(2, 2)
nargoutchk(1, 1)
training_velocity.time = 0;
validatevelocity(training_velocity)

% Extract data
data = permute(training_velocity.input, [1 3 2]);
teacher = training_velocity.sources;
n_sources = size(data, 1);

% Initialize ELM
elm.n_nodes = params.elm.n_nodes;
elm.f_internal = params.elm.f_internal;
elm.n_outputs = size(teacher, 2);
elm.n_inputs = size(data, 2);
elm.weights_internal = 2 * (0.5 - rand(elm.n_nodes, elm.n_inputs));
elm.bias_internal = 2 * (0.5 - rand(1, elm.n_nodes));

% Iterations to prevent memory issues
one_sample_size = elm.n_nodes * 8;
step_size = floor(params.memory.elm_iteration_limit / one_sample_size);
starts = 1:step_size:n_sources;

if step_size < 1.5 * elm.n_nodes
    warning('ELM:notEnoughSources',...
        'The initialisation step is too small for the number of nodes');
end

% Training loop
for start = starts
    indices = start:min(start+(step_size-1), n_sources);
    
    V = data(indices, :);
    T = teacher(indices, :);
    H = elm.f_internal((elm.weights_internal * V')' + elm.bias_internal);

    if start == starts(1)
        K = H'*H;
        B = (K\H') * T;
    else
        K = K + H'*H;
        B = B + (K\H') * (T - H*B);
    end
end
% Write output weights
elm.weights_output = B;

end
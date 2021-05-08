function [predictions, E] = predict_lcmv(velocity, templates, params)
%predict_knn
%
% Syntax: predictions = predict_lcmv(velocity, templates, params)
%
% Implements lcmv predictions.
narginchk(3, 3)
nargoutchk(1, 2)
validatevelocity(velocity)
validatevelocity(templates)

n_samples = size(velocity.input, 1);
n_templates = size(templates.input, 1);
n_times = size(velocity.input, 2);
predictions = zeros(n_samples, 3);

locations = templates.sources;
templates = permute(templates.input, [1 3 2]);
signals = permute(velocity.input, [3 2 1]);

E = zeros(n_samples, n_templates);

for sample_idx = 1:n_samples
    R = signals(:, :, sample_idx) * signals(:, :, sample_idx)';
    R = R ./ n_times; % does not change prediction
    E(sample_idx,:) = 1./sum(templates / R .* templates, 2);
    [~,i] = max(E(sample_idx,:));
    predictions(sample_idx,:) = locations(i,:);
end

velocity = apply_input_mode(...
    extract_signal_over_sensors(velocity, params),...
    params...
);

predictions = determine_orientation_phase(...
    predictions,...
    velocity,...
    params...
);

end
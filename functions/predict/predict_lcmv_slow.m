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
n_sensors = size(velocity.input, 3);
predictions = zeros(n_samples, 3);

locations = templates.sources;
templates = permute(templates.input, [3 2 1]);
signals = permute(velocity.input, [3 2 1]);

E = zeros(n_samples, n_templates);

for sample_idx = 1:n_samples
    R = zeros(n_sensors);
    s = signals(:, :, sample_idx);
    for time_idx = 1:n_times
        R = R + (s(:, time_idx) * s(:, time_idx)');
    end
    R = inv(R ./ n_times);
    
    e = zeros(1, n_templates);
    for template_idx = 1:n_templates
        d = templates(:, :, template_idx);
        e(template_idx) = 1./(d' * R * d);
    end
    E(sample_idx,:) = e;
    [~,i] = max(e);
    predictions(sample_idx,:) = locations(i,:);
end

end
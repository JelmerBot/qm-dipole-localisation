function [predictions, indices] = predict_knn(velocity, templates, params)
%predict_knn
%
% Syntax: [predictions, indices] = predict_knn(velocity, templates, params)
%
% Implements knn predictions.
narginchk(3, 3)
nargoutchk(1, 2)
velocity.time = 0; % force that time is removed
velocity.time = 0; % force that time is removed
validatevelocity(velocity)
validatevelocity(templates)

k = params.knn.k;
n_rows = size(velocity.input, 1);
indices = zeros(n_rows, k);
predictions = zeros(n_rows, 3);

locations = templates.sources;
templates = templates.input ./ max(templates.input, [], 3);
signals = velocity.input ./ max(velocity.input, [], 3);

parfor row = 1:n_rows
    difference = sqrt(sum((signals(row, :, :) - templates).^2, 3));
    [~, i] = mink(difference, k);
    predictions(row, :) = [...
        mean(locations(i, 1:2), 1),...
        compute_angle_average(locations(i, 3)')...
    ];
    indices(row,:) = i;
end

end
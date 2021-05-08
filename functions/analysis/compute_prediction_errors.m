function errors = compute_prediction_errors(predictions, teacher, params)
%compute_prediction_error
%
% Syntax: errors = compute_prediction_errors(predictions, teacher, params)
%
% Computes a normalised error value to optimize hyper parameters with
% the metric is unitless.

narginchk(3, 3)
nargoutchk(1, 1)
if any(size(predictions) ~= size(teacher))
    error('Predicitons and teacher need to be the same size');
end

d_x = predictions(:, 1) - teacher(:, 1);
d_y = predictions(:, 2) - teacher(:, 2);
d_azi = compute_angle_difference(teacher(:, 3), predictions(:, 3));

domain = params.domain;
range_x = diff(domain.x_range);
range_y = diff(domain.y_range);
range_azi = diff(domain.azimuth_range);

errors = sum(abs([d_x ./ range_x, d_y / range_y, d_azi / range_azi]), 2);

end

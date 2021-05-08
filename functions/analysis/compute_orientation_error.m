function errors = compute_orientation_error(predictions, teacher)
%compute_orientation_error
%
% Syntax: errors = compute_orientation_error(predictions, teacher)
%
% Computes the orientation error in radian.

narginchk(2, 2)
nargoutchk(1, 1)
if any(size(predictions) ~= size(teacher))
    error('Predicitons and teacher need to be the same size');
end

errors = abs(compute_angle_difference(teacher(:, 3), predictions(:, 3)));

end
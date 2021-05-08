function errors = compute_location_error(predictions, teacher)
%compute_location_error
%
% Syntax: errors = compute_location_error(predictions, teacher)
%
% Computes the location error in meter.

narginchk(2, 2)
nargoutchk(1, 1)
if any(size(predictions) ~= size(teacher))
    error('Predicitons and teacher need to be the same size');
end

errors = sqrt(sum((predictions(:, 1:2) - teacher(:, 1:2)).^2, 2));

end
function estimate = update_nonlin_estimate(estimate, direction, step_size, norm_limit, params)
%update_estimate
%
% Syntax:  estimate = update_nonlin_estimate(estimate, direction, step_size, norm_limit, params)
%
% Updates the estimate given a direction.
% - Clip the gradient to the given norm limit
% - Perform the update with the given (fixed) step size
% - Project the new estimate to the domain
narginchk(5, 5)
nargoutchk(1, 1)
validateattributes(estimate, {'double'}, {'nrows', 1, 'ncols', 3}, 'update_estimate', 'estimate', 1)
validateattributes(direction, {'double'}, {'nrows', 3, 'ncols', 1}, 'update_estimate', 'direction', 2)
validateattributes(step_size, {'double'}, {'scalar', 'positive', 'nonzero'}, 'update_estimate', 'step_size', 3)
validateattributes(norm_limit, {'double'}, {'scalar', 'positive', 'nonzero'}, 'update_estimate', 'norm_limit', 4) 

% Clip the norm of the direction
if norm(direction) > norm_limit
    direction = direction .* norm_limit / norm(direction);
end

% Perform step
estimate = estimate + step_size .* direction';

% Project to domain bounds
estimate(1) = clip(estimate(1), min(params.domain.x_range), max(params.domain.x_range));
estimate(2) = clip(estimate(2), min(params.domain.y_range), max(params.domain.y_range));
estimate(3) = mod(estimate(3), max(params.domain.azimuth_range));
    
end

function res = clip(v, v_min, v_max)
    res = min([max([v, v_min]), v_max]);
end
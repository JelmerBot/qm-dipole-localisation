function grad = compute_gradient_numerical(estimate, params, delta)
%compute_gradient_numerical
%
% Syntax:  grad = compute_gradient_numerical(estimate, params, delta)
%
% Numerically estimates the gradient of vx and vy over the x, y, and 
% azimuth dimensions. Delta is the step size in the estimation.
% 
% Output size is 2n x 3, where n is the number of sensors. The first
% n rows are for vx the last rows are for vy.
narginchk(3, 3)
nargoutchk(1, 1)
validateattributes(estimate, {'double'}, {'numel', 3}, 'compute_gradient_numerical', 'estimate', 1)

% Compute locations around the estimate
x_values = [estimate(1) - delta, estimate(1), estimate(1) + delta];
y_values = [estimate(2) - delta, estimate(2), estimate(2) + delta];
azi_values = [estimate(3) - delta, estimate(3), estimate(3) + delta];
[Y, X, A1] = ndgrid(y_values, x_values, azi_values);
grad_sources = [X(:) Y(:) A1(:)];
clear X Y A1

% Compute velocities around theta
velocity = compute_sensor_velocity(grad_sources, 0, params);
velocity = apply_input_mode(velocity, params);
x = reshape(abs(velocity.input), [3 3 3 size(velocity.input, 3)]);
[dx, dy, da] = gradient(x, delta, delta, delta, 1);
grad = [squeeze(dx(2,2,2,:)), squeeze(dy(2,2,2,:)), squeeze(da(2,2,2,:))];

end
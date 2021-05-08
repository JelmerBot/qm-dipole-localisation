function [predictions, estimates] = compute_nonlin_prediction(velocity, predict, max_iter)
%compute_nonlin_prediction
%
% Syntax: [predictions, estimates] = ...
%                 compute_nonlin_prediction(velocity, predict, max_iter)
%
% Computes a non-linear estimation of [x, y, azimuth] which fits best with
% the source velocity in signal. The velocities in signal should not have a
% time dimension.
%
% The predict argument should be a function handle to either the GN or 
% NR prediction function. Note, params should be captured by the handle
% It should have the form:
%   [prediction, estimates, iterations] = predict(M, estimates)

narginchk(3, 3)
nargoutchk(1, 2)
velocity.time = 0;              % Force that time dimension is removed!
validatevelocity(velocity)
validateattributes(predict, {'function_handle'}, {}, 'compute_nonlin_prediction', 'predict', 2)

% Only use the absolute signal
velocity.input = abs(velocity.input);

% Prepare to store the estimates
n_locations = size(velocity.xt, 1);
estimates = nan(n_locations, max_iter, 3);
predictions = nan(n_locations, 3);

% Loop over the source locations
input = velocity.input;
parfor row = 1:n_locations
   [predictions(row,:), estimates(row, :, :)] = ...
       predict(squeeze(input(row, :, :)), squeeze(estimates(row,:,:))); 
end

end
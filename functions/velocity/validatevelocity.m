function validatevelocity( velocity )
%validatevelocity

% This function is a validation of the data-structures used in the program
% Skip the checks when the program runs deployed. This saves time on the
% cluster
if isdeployed
    return
end

nargoutchk(0, 0);
narginchk(1, 1);

% For the parameter
validateattributes(velocity, {'struct'}, {'nonempty'}, mfilename, 'velocity', 1)
validateattributes(velocity.sources, {'numeric'}, {'real', '2d', 'size', [nan, 3]}, mfilename, 'velocity.sources', 1)
n_sources = size(velocity.sources, 1);
validateattributes(velocity.time, {'numeric'}, {'real', 'vector', 'increasing'}, mfilename, 'velocity.time', 1);

if ~isfield(velocity, 'power_xt')
    n_samples = length(velocity.time);
    validateattributes(velocity.xt, {'numeric'}, {'real', '3d', 'size', [n_sources, n_samples, nan]}, mfilename, 'velocity.xt', 1)
    validateattributes(velocity.yt, {'numeric'}, {'real', '3d', 'size', [n_sources, n_samples, nan]}, mfilename, 'velocity.yt', 1)
else
    n_samples = length(velocity.freqs);
    validateattributes(velocity.xt, {'numeric'}, {'3d', 'size', [n_sources, 2*n_samples, nan]}, mfilename, 'velocity.xt', 1)
    validateattributes(velocity.yt, {'numeric'}, {'3d', 'size', [n_sources, 2*n_samples, nan]}, mfilename, 'velocity.yt', 1)
    validateattributes(velocity.power_xt, {'numeric'}, {'real', '3d', 'size', [n_sources, nan, n_sensors]}, mfilename, 'velocity.power_xt', 1)
    validateattributes(velocity.power_yt, {'numeric'}, {'real', '3d', 'size', [n_sources, nan, n_sensors]}, mfilename, 'velocity.power_yt', 1)
end
    

end


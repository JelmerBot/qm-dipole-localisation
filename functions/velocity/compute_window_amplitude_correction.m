function factor = compute_window_amplitude_correction(params)
%compute_window_amplitude_correction
%
% Syntax: factor = compute_window_amplitude_correction(params)
%
% Computes the amplitude correction factor for the given window.
% See: 
% https://mathworks.com/matlabcentral/answers/372516-calculate-windowing-correction-factor
% for details
narginchk(1, 1)
nargoutchk(1, 1)
validateattributes(params.window, {'function_handle'}, {}, 'compute_window_amplitude_correction', 'window', 1)

w = params.window(params.sensors.sampling_rate);
factor = length(w) / sum(w);
end


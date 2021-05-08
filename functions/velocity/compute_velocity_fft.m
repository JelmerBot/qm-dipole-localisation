function velocity = compute_velocity_fft(velocity, params)
%compute_velocity_fft
%
% Syntax: velocity = compute_velocity_fft(velocity, params)
%
% Computes the DFT power.
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(velocity)
validateattributes(params.window, {'function_handle'}, {}, 'compute_velocity_fft', 'window', 2)

Fs = params.sensors.sampling_rate;
Fi = params.source.frequency;
M = size(velocity.xt, 2);
window_correction = compute_window_amplitude_correction(params);

velocity.freqs = Fs * (0:floor(M/2)) / M;
[~, velocity.source_freq_idx] = min(abs(velocity.freqs - Fi));

% Compute spectrum
velocity.xt = fft(params.window(M)' .* velocity.xt, [], 2);
velocity.yt = fft(params.window(M)' .* velocity.yt, [], 2);

power = abs(velocity.xt / M);
power = power(:, 1:floor(M/2)+1, :);
power(:, 2:end-1,:) = 2 * power(:, 2:end-1,:);
velocity.power_xt = window_correction .* power;

power = abs(velocity.yt / M);
power = power(:, 1:floor(M/2)+1, :);
power(:, 2:end-1,:) = 2 * power(:, 2:end-1,:);
velocity.power_yt = window_correction .* power;

end



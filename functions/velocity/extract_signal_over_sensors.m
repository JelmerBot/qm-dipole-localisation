function velocity = extract_signal_over_sensors(velocity, params)
%extract_signal_over_sensors
%
% Syntax: velocity = extract_signal_over_sensors(velocity, params)
%
% Extracts a time averaged velocity signal over the sensors.

narginchk(2, 3)
nargoutchk(1, 1)
validatevelocity(velocity)
validateattributes(params.window, {'function_handle'}, {}, 'compute_velocity_fft', 'window', 2)

% Compute fft
velocity = compute_velocity_fft(velocity, params);

% Extract source frequency
i_source = velocity.source_freq_idx;

phase = angle(velocity.xt(:, i_source, :));
phase = wrapToPi(phase);
amplitude = velocity.power_xt(:, i_source, :);
velocity.xt = sign(phase) .* amplitude; 

phase = angle(velocity.yt(:, i_source, :));
phase = wrapToPi(phase) / pi;
amplitude = velocity.power_yt(:, i_source, :);
velocity.yt = sign(phase) .* amplitude; 

velocity.time = 0;

velocity = rmfield(velocity, 'power_xt');
velocity = rmfield(velocity, 'power_yt');
velocity = rmfield(velocity, 'freqs');
velocity = rmfield(velocity, 'source_freq_idx');

end


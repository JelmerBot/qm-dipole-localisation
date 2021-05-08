function time = compute_time(params)
%compute_time
%
% Syntax: time = compute_time(params)
%
% How the time is sampled is really important for the accuracy of the pase
% when computed with an FFT. So don't change this. If you do, make sure to
% check whether the pase of the extract_signal_over_sensors is still
% accurate!
narginchk(1, 1);
nargoutchk(1, 1);

sensors = params.sensors;
domain = params.domain;
total_length = domain.measurement_duration * sensors.sampling_rate;
period = 1 / sensors.sampling_rate;
time = (1:total_length) * period;

end
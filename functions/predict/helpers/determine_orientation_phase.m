function predictions = determine_orientation_phase(predictions, velocity, params)
%determine_orientation_phase
%
% Syntax: predictions = determine_orientation_phase(predictions, velocity, params)
%
% Some predictors cannot distinghuis the phase of a dipole (GN, NR, LCMV).
% So match both the estimated orientation and the estimate + pi to the velocity
% measurements to figure out which one is best.


estimated_velocity = apply_input_mode(...
    compute_sensor_velocity(predictions, 0, params), params...
);

estimated_signal = estimated_velocity.input;
flipped_signal = -estimated_velocity.input;
error_estimated_signal = sqrt(sum((estimated_signal - velocity.input).^2, 3));
error_flipped_signal = sqrt(sum((flipped_signal - velocity.input).^2, 3));

% Update the predictions which match better when flipped.
require_update = error_estimated_signal > error_flipped_signal;
predictions(:, 3) = mod(predictions(:, 3) + require_update * pi, 2*pi);

end


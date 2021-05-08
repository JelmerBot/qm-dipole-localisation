function velocity = simulate_lsq_velocity(source, params)
%simulate_lsq_velocity
%
% Syntax: velocity = simulate_lsq_velocity(source, params)
%
% source(1) = x, source(2) = y, source(3) = azi

velocity = compute_sensor_velocity(source, 0, params);
velocity = apply_input_mode(velocity, params);
velocity = permute(velocity.input, [1 3 2]);

end
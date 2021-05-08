function velocity = apply_sensor_noise(velocity, params)
%apply_sensor_noise
%
% Syntax: velocity = apply_sensor_noise(velocity, params)

narginchk(2, 2);
nargoutchk(1, 1);

sensors = params.sensors;
velocity.xt = velocity.xt + sensors.sensing_noise .* randn(size(velocity.xt));
velocity.yt = velocity.yt + sensors.sensing_noise .* randn(size(velocity.yt));

end
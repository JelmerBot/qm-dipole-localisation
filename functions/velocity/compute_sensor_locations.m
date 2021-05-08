function locations = compute_sensor_locations(sensors)
%compute_sensor_locations
%
% Syntax: locations = compute_sensor_locations(sensors)

narginchk(1, 1);
nargoutchk(1, 1);

locations.x = linspace(min(sensors.x_range), max(sensors.x_range), sensors.n_sensors);
locations.y = repmat(sensors.y_value, 1, sensors.n_sensors);

end
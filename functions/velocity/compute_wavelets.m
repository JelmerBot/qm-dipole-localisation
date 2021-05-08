function wavelets = compute_wavelets(sources, params)
%compute_wavelets
%
% Syntax: wavelets = compute_wavelets(sources, params)
narginchk(2, 2)
nargoutchk(1, 1)

source = params.source;
sensors = params.sensors;
domain = params.domain;

% Positions
step = diff(sensors.locations.x(1:2));
x1 = fliplr((sensors.locations.x(1)-step):-step:min(domain.x_range));
x3 = (sensors.locations.x(end)+step):step:max(domain.x_range);
x = [x1 sensors.locations.x x3];
y = sensors.y_value;
b = sources(:,1);
d = sources(:,2);
sigma = (x - b) ./ (y - d);

% Velocity
w = 2 .* pi .* source.frequency .* source.amplitude .*...
    cos(2 .* pi .* source.frequency);
scaling = (source.radius.^3 .* w) ./ (2 .* abs(y - d).^3);

% Wavelets
denom = (sigma.^2 + 1).^(5/2);
even = scaling .* (2.*sigma.^2 - 1) ./ denom;
odd =  scaling .* (3.*sigma) ./ denom;
navelet = scaling .* (2 - sigma.^2) ./ denom;

% Output (normalised wavelets)
wavelets.x = x;
sensor_mask = [false(1, length(x1)), true(1, sensors.n_sensors), false(1, length(x3))];
wavelets.sources = sources;
wavelets.even = even(:, sensor_mask) ./ max(abs(even), [], 2);
wavelets.odd = odd(:, sensor_mask) ./ max(abs(odd), [], 2);
wavelets.navelet = navelet(:, sensor_mask) ./ max(abs(navelet), [], 2);

end

function wavelets = compute_wavelets_raw(sources, params)
%compute_wavelets
%
% Syntax: wavelets = compute_wavelets_raw(sources, params)
narginchk(2, 2)
nargoutchk(1, 1)

source = params.source;
sensors = params.sensors;

% Positions
x = sensors.locations.x;
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

% Output (normalised wavelets)
wavelets.x = x;
wavelets.sources = sources;
wavelets.even = scaling .* (2.*sigma.^2 - 1) ./ denom;
wavelets.odd = scaling .* (3.*sigma) ./ denom;
wavelets.navelet = scaling .* (2 - sigma.^2) ./ denom;

end

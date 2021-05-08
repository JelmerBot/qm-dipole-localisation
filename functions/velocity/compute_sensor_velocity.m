function velocity = compute_sensor_velocity(sources, time, params)        
%compute_sensor_velocity
%
% Syntax: velocity = compute_sensor_velocity(sources, time, params)       
% 
% Implements potential flow formulas. See the potential flow chapter for
% the formulas.

narginchk(3, 3);
nargoutchk(1, 1);

% Easy access variables
x = sources(:,1);
y = sources(:,2);
azi = sources(:,3);
source = params.source;
sensors = params.sensors;

% Relative sensor position from source
r_xt = permute(repmat(sensors.locations.x - x(:), 1, 1, length(time)), [1 3 2]);
r_yt = permute(repmat(sensors.locations.y - y(:), 1, 1, length(time)), [1 3 2]);

% Dispacement over time
r_t = source.amplitude .* sin(2 .* pi .* source.frequency .* time);
% Relative sensor positions from source over time
% use a  '-' here because it is 'r = sensor - source'
% so a positive displacement of the source reduces r
r_xt = r_xt - cos(azi(:)) .* r_t;
r_yt = r_yt - sin(azi(:)) .* r_t;
r_mt = sqrt(r_xt.^2 + r_yt.^2);

% Velocity over time
w_t = 2 .* pi .* source.frequency .* source.amplitude .*...
      cos(2 .* pi .* source.frequency .* time);
w_xt = cos(azi(:)) .* w_t;
w_yt = sin(azi(:)) .* w_t;

% Compute sensor velocity
dot = w_xt .* r_xt + w_yt .* r_yt;

velocity.xt = ...
    source.radius.^3 ./ (2 .* r_mt.^3) .* ...
    ( -w_xt + 3 .* r_xt .* dot ./ r_mt.^2);
velocity.yt = ...
    source.radius.^3 ./ (2 .* r_mt.^3) .* ...
    ( -w_yt + 3 .* r_yt .* dot ./ r_mt.^2);

velocity.time = time;
velocity.sources = [x(:), y(:), azi(:)];

% Check if the contract holds
validatevelocity(velocity);

end

function predictions = predict_random(velocity, params)
%predict_random
%
% Syntax: predictions = predict_random(velocity, params)
%
% Implements uniform random predictions over the domain.
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(velocity)

domain = params.domain;
n_sources = size(velocity.xt, 1);
x_values = uniform_interval([n_sources, 1], min(domain.x_range), max(domain.x_range));
y_values = uniform_interval([n_sources, 1], min(domain.y_range), max(domain.y_range));
azi_values = uniform_interval([n_sources, 1], min(domain.azimuth_range), max(domain.azimuth_range));

predictions = [x_values, y_values, azi_values];

end

function values = uniform_interval(sz, min, max)
    values = min + (max - min) .* rand(sz);
end
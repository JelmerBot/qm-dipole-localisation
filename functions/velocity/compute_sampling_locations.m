function source_info = compute_sampling_locations(params)
%compute_sampling_locations - Generates sampling grid
%
% Syntax: source_info = compute_sampling_locations(params)
    
narginchk(1, 1)
nargoutchk(1, 1)

domain = params.domain;

% Find the domain sizes
l_x = diff(domain.x_range);
if l_x == 0; l_x = []; end
l_y = diff(domain.y_range);
if l_y == 0; l_y = []; end
l_az = diff(domain.azimuth_range) / (2*pi); 
if l_az == 0; l_az = []; end

scale = 1000;

% Sample sources in the domain
sz = scale .* [l_x l_y l_az];
if ~isempty(sz)
    source_info = poissonDisc(sz, scale .* domain.min_source_distance);
    source_info = source_info ./ scale;
else
    error('Error generating locations')
end

% Some dimensions may be empty check them all
if isempty(l_x)
    source_info = [repmat(domain.x_range, size(source_info, 1), 1) source_info];
else
    source_info(:, 1) = source_info(:, 1) + min(domain.x_range);
end
if isempty(l_y)
    source_info = [source_info(:, 1), repmat(domain.y_range, size(source_info, 1), 1), source_info(:, 2:end)];
else
    source_info(:, 2) = source_info(:, 2) + min(domain.y_range);
end
if isempty(l_az)
    source_info = [source_info(:, 1:2), repmat(domain.azimuth_range, size(source_info, 1), 1), source_info(:,3:end)];
else
    source_info(:, 3) = source_info(:, 3) * 2*pi + min(domain.azimuth_range);
end
end
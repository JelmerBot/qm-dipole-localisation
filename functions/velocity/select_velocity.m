function [velocity, idx] = select_velocity(velocity, pos, phi)
%select_velocity - selects from velocity samples
%
% Syntax: velocity = select_velocity(velocity, pos, phi)
narginchk(3, 3)
nargoutchk(1, 2)
validatevelocity(velocity);

idx = 1:size(velocity.xt, 1);

if ~isempty(pos) && isempty(phi)
    validateattributes(pos, {'numeric'}, {'real', 'row', 'ncols', 2}, mfilename, 'pos', 2)
    dist = sum((pos - velocity.sources(:, 1:2)).^2, 2);
    idx = find(dist == min(dist));
elseif isempty(pos) && ~isempty(phi)
    validateattributes(phi, {'numeric'}, {'real', 'scalar'}, mfilename, 'phi', 3)
    dist = sum((phi - velocity.sources(:, 3)).^2, 2);
    idx = find(dist == min(dist));
elseif ~isempty(pos) && ~isempty(phi)
    validateattributes(pos, {'numeric'}, {'real', 'row', 'ncols', 2}, mfilename, 'pos', 2)
    validateattributes(phi, {'numeric'}, {'real', 'scalar'}, mfilename, 'phi', 3)
    dist = sum(([pos phi] - velocity.sources).^2, 2);
    idx = find(dist == min(dist));
end

velocity = select_velocity_indices(velocity, idx);

end
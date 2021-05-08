function velocity = select_velocity_indices(velocity, indices)
%select_velocity_indices
%
% Syntax: velocity = select_velocity_indices(velocity, indices)
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(velocity);

velocity.xt = velocity.xt(indices, :, :);
velocity.yt = velocity.yt(indices, :, :);
velocity.sources = velocity.sources(indices,:);
if isfield(velocity, 'input')
    velocity.input = velocity.input(indices,  :, :);
end
if isfield(velocity, 'power_xt')
    velocity.power_xt = velocity.power_xt(indices,:,:);
end
if isfield(velocity, 'power_yt')
    velocity.power_yt = velocity.power_yt(indices,:,:);
end

end
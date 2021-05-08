function velocity = normalise_components(velocity, params)
%normalise_components
%
% Syntax: velocity = normalise_components(velocity, params)
narginchk(1, 2)
nargoutchk(1, 1)
validatevelocity(velocity)

if strcmp(params.sensors.input_mode, 'x|y')
    velocity.xt = velocity.xt ./ max(abs(velocity.xt(:, :, 1:2:end)), [], 3);
    velocity.yt = velocity.yt ./ max(abs(velocity.yt(:, :, 2:2:end)), [], 3);
else
    velocity.xt = velocity.xt ./ max(abs(velocity.xt), [], 3);
    velocity.yt = velocity.yt ./ max(abs(velocity.yt), [], 3);
end

end
function quadrature = simulate_quadrature(source, params)
%simulate_quadrature
%
% Syntax: quadrature = simulate_quadrature(x, params)
%
% source(1) = x, source(2) = y, source(3) = azi

velocity = compute_sensor_velocity(source, 0, params);

quadrature = sqrt(velocity.xt.^2 + 0.5 .* velocity.yt.^2);
quadrature = permute(quadrature, [1 3 2]);

end
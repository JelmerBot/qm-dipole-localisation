function med = compute_angle_median(angles)
%compute_angle_median
%
% Syntax: med = compute_angle_median(angles)

narginchk(1, 1)
nargoutchk(1, 1)

% Angles between 0 and pi
angles = mod(angles, 2*pi);
med = mod(atan2(median(sin(angles), 2), median(cos(angles), 2)), 2*pi);

end
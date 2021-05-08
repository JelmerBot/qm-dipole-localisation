function average = compute_angle_average(angles)
%compute_angle_average
%
% Syntax: average = compute_angle_average(angles)
%
% Ignores the phase of the signal so averages range from 0 to pi
% Two averages are possible: a horizontal and vertical average.
% The average for which the given angles are closest together is returned
% 
% See https://en.wikipedia.org/wiki/Mean_of_circular_quantities for
% detail about computing an average angle
narginchk(1, 1)
nargoutchk(1, 1)

% Angles between 0 and 2*pi
angles = mod(angles, 2*pi);
average = mod(atan2(mean(sin(angles), 2), mean(cos(angles), 2)), 2*pi);

end

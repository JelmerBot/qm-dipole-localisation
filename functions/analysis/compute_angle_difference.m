function diff = compute_angle_difference(angle_1, angle_2)
%compute_angle_difference
%
% Syntax: diff = compute_angle_difference(angle_1, angle_2)
%
% Range from -pi to pi
narginchk(2, 2)
nargoutchk(1, 1)

% Rounds small values down to zero. 
% pol=sign(sin(angle_2-angle_1));
% diff=pol .* acos(1 - ((cos(angle_1) - cos(angle_2)).^2 +...
%                       (sin(angle_1) - sin(angle_2)).^2) / 2);

diff = atan2(sin(angle_2 - angle_1), cos(angle_2 - angle_1));
                  
end

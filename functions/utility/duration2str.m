function [ str ] = duration2str( duration_s )
%duration2str
%
% Syntax: str = duration2str(duration_s)
%
% The input should be in seconds!

narginchk(1, 1)
nargoutchk(0, 1)
if length(duration_s) > 1
    error('input should be a single duration in seconds');
end

% seconds
if duration_s < 60
    str = [num2str(duration_s), 's'];
    return;
end

% mins
duration_m = floor(duration_s / 60);
duration_s = mod(round(duration_s), 60);
if duration_m < 60
    str = [num2str(duration_m), 'm ',...
           num2str(duration_s), 's'];
    return;
end

% hours
duration_h = floor(duration_m / 60);
duration_m = mod(duration_m, 60);
if duration_h < 24
    str = [num2str(duration_h), 'h ',...
           num2str(duration_m), 'm ',...
           num2str(duration_s), 's'];
    return;
end

% days
duration_d = floor(duration_h / 24);
duration_h = mod(duration_h, 24);

str = [num2str(duration_d), 'd ',...
       num2str(duration_h), 'h ',...
       num2str(duration_m), 'm ',...
       num2str(duration_s), 's'];

end


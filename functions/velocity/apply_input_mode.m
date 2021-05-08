function velocity = apply_input_mode(velocity, params)
%apply_input_mode
%
% input = apply_input_mode(velocity, params)
%
% Dim1 = sources
% Dim2 = time
% Dim3 = sensors + component
narginchk(2, 2)
nargoutchk(1, 1)
validatevelocity(velocity)

% All sensors measure both x and y
if strcmp(params.sensors.input_mode, 'x+y')
    % nsensors * x nsensors * y
    velocity.input = cat(3, velocity.xt, velocity.yt);
    return
end
   
% All sensors measure x
if strcmp(params.sensors.input_mode, 'x')
    % nsensors * x
    velocity.input = velocity.xt;
    return
end

% All sensors measure y
if strcmp(params.sensors.input_mode, 'y')
    % nsensors * y
    velocity.input = velocity.yt;
    return
end

% Subsequent sensors measure x then y
if strcmp(params.sensors.input_mode, 'x|y')
    i1 = 1:2:(params.sensors.n_sensors);
    i2 = 2:2:params.sensors.n_sensors;
    velocity.input = zeros(size(velocity.xt));
    velocity.input(:,:,i1) = velocity.xt(:,:,i1);
    velocity.input(:,:,i2) = velocity.yt(:,:,i2);
    return
end

error('unknown input mode');
end

function [vx, vy] = cwt_input_mode_velocity(velocity, params)

% All sensors measure both x and y
if strcmp(params.sensors.input_mode, 'x+y')
    vx = velocity.xt;
    vy = velocity.yt;
    return
end

% All sensors measure x
if strcmp(params.sensors.input_mode, 'x')
    vx = velocity.xt;
    vy = zeros(size(vx));
    return
end

% All sensors measure x
if strcmp(params.sensors.input_mode, 'y')
    vy = velocity.yt;
    vx = zeros(size(vy));
    return
end

% Subsequent sensors measure x then y
if strcmp(params.sensors.input_mode, 'x|y')
    i1 = 1:2:(params.sensors.n_sensors);
    i2 = 2:2:params.sensors.n_sensors;
    
    vx = velocity.xt;
    vy = velocity.yt;
    vx(:,:,i2) = 0; % disable alternating sensors
    vy(:,:,i1) = 0; % disable alternating sensors
    return
end

end
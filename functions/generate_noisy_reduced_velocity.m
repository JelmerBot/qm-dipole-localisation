function output = generate_noisy_reduced_velocity(sources, params)
%generate_noisy_reduced_velocity
%
%   Syntax output = generate_noisy_reduced_velocity(sources, params)
%
% Simulates the velocity at the sensors using:
%   extract_signal_over_sensors(...
%       apply_sensor_noise(...
%           compute_sensor_velocity(sources,...)...
%       )...
%   )
% The output is computed in iterations to maxize memory usage without
% overloading the system (either the cluster at 64G or a desktop with 16G).

% Determine iteration steps
memory = params.memory;
domain = params.domain;
sensors = params.sensors;
n_time_samples = domain.measurement_duration * sensors.sampling_rate;
one_sample = (3 * sensors.n_sensors * n_time_samples) * 8;
stepsize = floor(memory.sim_iteration_limit / one_sample);
n_sources = size(sources, 1);
starts = 1:stepsize:n_sources;
n_iterations = length(starts);

if params.print_messages
    disp('Start velocity generation')
    disp(['    using ', num2str(n_iterations), ' iterations for this setting']);
end

% Preallocate for the results
time = compute_time(params);
output.xt = zeros(n_sources, 1, sensors.n_sensors);
output.yt = zeros(n_sources, 1, sensors.n_sensors);
output.sources = sources;
output.time = 0;

% Run the iterations
for idx = starts
    tic;
    indices = idx:min(idx+(stepsize-1), n_sources);
    
    velocity = compute_sensor_velocity(sources(indices,:), time, params);
    velocity = apply_sensor_noise(velocity, params);
    velocity = extract_signal_over_sensors(velocity, params);
    output.xt(indices, :, :) = velocity.xt;
    output.yt(indices, :, :) = velocity.yt;
    
    if params.print_messages && idx == 1
        t = toc;
        disp(['    iteration in ', duration2str(t)]);
        disp(['        estimated time: ', duration2str(n_iterations * t)])
    end
end

% Check postconditions
validatevelocity(output)

end
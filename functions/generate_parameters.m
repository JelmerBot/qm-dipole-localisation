function params = generate_parameters()
%generate_parameters - Generates default parameters
%
% Syntax: params = generate_parameters()

narginchk(0, 0)
nargoutchk(1, 1)

% Control the program's output
params.write_output = true;
params.print_messages = true;
params.show_plots = ~isdeployed && false;

% Describe the sensors. n_sensors and sensing_noise may be 
% overwritten by a parameftersweep. See the sweep parameters.
sensors.n_sensors = 8;
sensors.x_range = [-0.2 0.2];     % x - locations (in m)
sensors.y_value = 0;              % y - location (in m)
sensors.sensing_noise = 1e-5;     % in m/s flow added to the sensors (Gaussian)
sensors.sampling_rate = 2048;     % in Hz
% x+y   All sensors measure x and y
% x     All sensors measure x
% y     All sensors measure y
% x|y   Subsequent sensors measure x then y
sensors.input_mode = 'x+y';       
sensors.locations = compute_sensor_locations(sensors);
params.sensors = sensors;

% Describe the source. A dipole is a vibrating sphere. 
source.amplitude = 2E-3;       % in m (is mm) centre-peak amplitude
source.frequency = 45;         % in Hz
source.radius = 1E-2;          % in m (is cm)
params.source = source;

% The 2D location domain and 1D motion domain. Together they form a 3D
% space of possible sources.
domain.x_range = [-0.5 0.5];      % along the sensor array (m)
domain.y_range = [0 0.5] +...     % perpendicular to the sensor array (m)
                 source.radius + 0.015; % at least 1.5cm between source and sensors
domain.azimuth_range = [0 2*pi];  % motion vector angle in the x-y plane
% Poisson disc sampling is used to generate source locations and motions
% This setting determines how far away sources should be from each other 
% in the 3D domain described above. The unit is m*rad (?). If two points have
% the same motion vector, they are at least x cm apart in location. If they
% have the same location they are at least 2*pi*x rad apart in motion
% direction.
domain.min_source_distance = 0.015; % value gets replaced by the experimental condition
domain.measurement_duration = 1;    % measurement duration in seconds
params.domain = domain;

% Configure training, validation, and testing
experiment.test_source_distance = 0.01; 
experiment.train_source_distances = [0.01 0.03 0.05 0.09];
experiment.input_modes = {'x+y', 'x', 'y', 'x|y'};
experiment.n_fold = 5;
params.experiment = experiment;

% Hyper-parameters for the elm.
elm.n_nodes = 800;
elm.f_internal = @relu; % tanh or relu
elm.validate.evaluations = 101;
params.elm = elm;

% Hyper-parameters for the cwt.
cwt.c_x = 25/45;
cwt.c_y = 45/52;
cwt.threshold_min = 0.1;
cwt.threshold_max = 0.9;
cwt.validate.evaluations = 30;
cwt.validate.threshold_min = [0, 1];
cwt.validate.threshold_max = [0, 1];
cwt.validate.c_x = [0.5 1];
cwt.validate.c_y = [0.3 0.9];
params.cwt = cwt;

% Hyper-parameters for the knn method.
knn.k = 2;
knn.validate.k = 1:20; % values for grid-search
params.knn = knn;

% Hyper-parameters for the gn method. See predict_gn
% for what each parameter value does.
gn.starting_estimate = [0 0.05 pi]; % x y azimuth
gn.max_iterations = 100;
gn.gradient_step = 1e-3;
gn.step_tolerance = 1e-3;
gn.step_size = 1;
gn.norm_limit = inf;
gn.validate.starting_y = domain.y_range;
gn.validate.evaluations = 21;
params.gn = gn;

% Hyper-parameters for the nr method. See help predict_nr
% for what each parameter value does.
nr.starting_estimate = [0 0.05 pi]; % x y azimuth
nr.max_iterations = 100;
nr.gradient_step = 1e-3;
nr.step_tolerance = 1e-3;
nr.step_size = 1;
nr.norm_limit = 0.2;
% nr.regularization = 1e-15;
nr.validate.starting_y = domain.y_range;
% nr.validate.step_size = [0.1 1];
nr.validate.norm_limit = [0.1 1];
% nr.validate.regularization = [eps 1e-10];
nr.validate.evaluations = 30;
params.nr = nr;

% Hyper-parameters for the quadrature method. 
qm.refine_repeats = 4;
qm.refine_iterations = 10;     
qm.gradient_step = 1e-3;       % Same as GN
qm.step_tolerance = 1e-3;      % Same as GN
qm.function_tolerance = eps;   % Almost disable this
qm.optimality_tolerance = eps; % Almost disable this
params.qm = qm;

% Hyper-paramters for the lsq method
lsq.starting_estimate = [0 0.05 pi];
lsq.max_iterations = 100;       % Same as GN
lsq.gradient_step = 1e-3;       % Same as GN
lsq.step_tolerance = 1e-3;      % Same as GN
lsq.function_tolerance = eps;   % Almost disable this
lsq.optimality_tolerance = eps; % Almost disable this
params.lsq = lsq;

% Hyper-parameters for the mlp method. 
% Validated parameters
mlp.layer_sizes = [512 512 512]; % validate
mlp.learn_rate = 0.001;          % validate
mlp.verbose = true;              % Toggle training output
% Fixed parameters
mlp.max_epochs = 500;                        % Hand validated
mlp.validation_patience = 5;                 % MATLAB default
mlp.training_levels = [0.2 0.4 0.6 Inf];     % Hand validated
mlp.input_normalisation = true;              % Problem statement
mlp.input_transform = 'none';                % Problem statement
mlp.initializer = @glorot;                   % MATLAB default (and literature)
mlp.use_dropout = false;                     % Not used (worse prelim results)
mlp.input_dropout = 0.2;                     % Based on literature
mlp.hidden_dropout = 0.5;                    % Based on literature
mlp.use_batch_normalisation = false;         % Not used (worse prelim results)
mlp.solver = 'adam';                         % Based on literature
mlp.error_metric = 'mean_absolute_error';    % Possible values: 'mean_squared_error', 'mean_absolute_error', 'median_absolute_error'
mlp.execution_environment = 'cpu';           % For the cluster
mlp.minibatch_size = 2048;                   % For the cluster
mlp.gradient_threshold = Inf;                % Not used
mlp.learn_decay_schedule = 10;               % MATLAB default
mlp.learn_decay_factor = 0.1;                % MATLAB default
% Configure validation
mlp.validate.n_layers = [1 4];
mlp.validate.layer_size = [16 1024];
mlp.validate.learn_rate = [1e-4 1e-1];
mlp.validate.evaluations = 30;
params.mlp = mlp;

% The implementation is memory bound. These values limit how much memory
% is used per iteration at several stages in the program. The values are
% tuned for a PC with 16G RAM when the code is run in matlab. When the code
% is deployed (compiled) the values are tuned for a high performance computer 
% with (64G RAM).
if isdeployed % cluster
    % The simulation computes the velocity and LCMV prediction of
    % multiple sources in one iteration. This value determines how many
    % sources are used in one iteration. Especially generating the
    % velocity over time and computing the corresponding pattern over the
    % sensors takes much memory. This parameter determines the size of
    % velocity.xt and velocity.yt;
    memory.sim_iteration_limit = 7.5e9;
    % The OS-ELM memory usage is bounded by a matrix K (or P) which is 
    % NxN where N is the number of nodes. The upperbound of N in the 
    % validation sweep depends on the number on unique sources. Generally,
    % the number of unqiue sources is way to high. So this parameter is
    % used as an additional upperbound. 
    memory.elm_nodes_limit = 4e9;
    % The second OS-ELM parameter determines how large matrix H can be. 
    % Based on that the number of rows is computed. If there are too few 
    % samples in an iteration H becomes singular and the inversion 
    % inacurate. If that is the case reduce elm_nodes_limit and increase
    % elm_iteration_limit. Make sure that elm_nodes_limit is large enough
    % that the validation sweep find the best number of nodes.
    memory.elm_iteration_limit = 17e9;
    memory.mem = 62;
else          % desktop
    memory.sim_iteration_limit = 1.9e9;
    memory.elm_nodes_limit = 4e9;
    memory.elm_iteration_limit = 1.5e9;
    memory.mem = 16;
end
c = parcluster('local');
memory.max_workers = c.NumWorkers;
params.memory = memory;

params.output_folder = make_sweep_name();

end

function name = make_sweep_name()
    if isdeployed
        base_folder = '/data/';
    else
        base_folder = './data/';
    end
    name = fullfile(base_folder, date);
end

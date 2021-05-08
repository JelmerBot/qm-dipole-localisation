%% Introduction
% Testing the CWT orientation prediction. 

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Analytical value of c_x and c_y

syms rho

denom = (rho.^2 + 1).^(5/2);

even = (2*rho.^2 - 1) ./ denom;
odd = (3*rho) ./ denom;
navelet = (2 - rho.^2) ./ denom;

int_even = int(even.^2, rho, [-Inf, Inf]);
int_odd = int(odd.^2, rho, [-Inf, Inf]);
int_nav = int(navelet.^2, rho, [-Inf, Inf]);

c_x = int_even ./ int_odd
c_y = int_odd ./ int_nav

%% Numerical value of c_x and c_y
% Using A43 and A48, we can compute c_x and c_y using sources with known
% locations and orientations

% Setup parameters
params = generate_parameters();
params.sensors.n_sensors = 64;
params.sensors.locations = compute_sensor_locations(params.sensors);
load('./config/parameters/window.mat');
params.window = window;

% Compute velocity and wavelets
b = linspace(-0.5, 0.5, 20);
d = [0.1 0.2 0.3 0.4];
phis = linspace(0, 2*pi, 20); 
[B, D, P] = meshgrid(b, d, phis);
sources = [B(:), D(:), P(:)];

velocity = generate_noisy_reduced_velocity(sources, params);
[vx, vy] = cwt_input_mode_velocity(velocity, params);
vx = permute(vx, [1 3 2]);
vy = permute(vy, [1 3 2]);

wavelets = compute_wavelets(sources(:, 1:2), params);
even = wavelets.even;
odd = wavelets.odd;
navelet = wavelets.navelet;

% Compute cwt coefficients for the known source location
scaling = 1 ./ sqrt(abs(params.sensors.y_value - wavelets.sources(:,2)));
x_even = scaling .* sum(vx .* even, 2);
x_odd =  scaling .* sum(vx .* odd, 2); 
y_nav =  scaling .* sum(vy .* navelet, 2);
y_odd =  scaling .* sum(vy .* odd, 2);

% Compute correction coefficients and predictions
c_x = tan(sources(:, 3)) ./ (x_odd ./ x_even);
c_y = tan(sources(:, 3)) ./ (y_nav ./ y_odd);

output = table(...
    sources(:, 1), sources(:, 2), sources(:, 3),...
    c_x, c_y, x_odd, x_even, y_odd, y_nav, ...
    'VariableNames', {'b', 'd', 'phi', 'c_x', 'c_y', 'x_odd', 'x_even', 'y_odd', 'y_nav'});

clear scaling x_even x_odd y_nav y_odd c_x c_y

%% Analyse numerical c_x and c_y values

[g, b, d] = findgroups(output.b, output.d);

% Median c_x over source orientation. Because source orientation is not
% known at the time of orientation prediction.
med_c_x = splitapply(@median, output.c_x, g);
mean_c_x = splitapply(@mean, output.c_x, g);
% Values of c_x for locations towards b = 0 at two distances
far_lim_x = median(med_c_x(d == 0.4 & b > -0.1 & b < 0.1));
close_lim_x = median(med_c_x(d == 0.1 & b > -0.1 & b < 0.1));

figure(3)
hold off
plot3(b, d, med_c_x, '.')
hold on
plot3(b, d, mean_c_x, '.')
plot3([-0.5 0.5], [0 0], [far_lim_x, far_lim_x], '--')
plot3([-0.5 0.5], [0 0], [close_lim_x, close_lim_x], '--')
zlim([-2 2])
grid on

% Median c_y over source orientation. Because source orientation is not
% known at the time of orientation prediction.
med_c_y = splitapply(@median, output.c_y, g);
mean_c_y = splitapply(@mean, output.c_y, g);
% Values of c_y for locations towards b = 0 at two distances
far_lim_y = median(med_c_y(d == 0.4 & b > -0.1 & b < 0.1));
close_lim_y = median(med_c_y(d == 0.1 & b > -0.1 & b < 0.1));

figure(4)
hold off
plot3(b, d, med_c_y, '.')
hold on
plot3(b, d, mean_c_y, '.')
plot3([-0.5 0.5], [0 0], [far_lim_y, far_lim_y], '--')
plot3([-0.5 0.5], [0 0], [close_lim_y, close_lim_y], '--')
zlim([-2 2])
grid on

%% Plot estimated orientation error for different locations

d = unique(output.d);
distance_str = arrayfun(@(idx) ['d=' num2str(d(idx)) ' m'], 1:length(d), 'UniformOutput', false);
c_xs = [1, 27/45, 25/45, 23/45];
c_x_label = {'1', '27/45', '25/45', '23/45'};
c_ys = [1, 45/123, 45/87, 45/52];
c_y_label = {'1', '45/123', '45/87', '45/52'};
linespec = {'-', '--', '-.', ':'};

figure(1)
colors = get(gca, 'ColorOrder');
n_rows = 2;
n_cols = length(d);


for col = 1:n_cols
    
    subplot(n_rows, n_cols, col);
    hold off
%     xlabel('$b$ (m)')
%     ylabel('RMSE of $\varphi$ (rad)')
%     ylim([0 1])
    subplot(n_rows, n_cols, col+n_cols)
    hold off
%     xlabel('$b$ (m)')
%     ylabel('RMSE of $\varphi$ (rad)')
%     ylim([0 1])
    
    values = output(output.d == d(col), :);
        
    for c = 1:length(c_xs)
        [b, e_x, e_y] = phi_error(values, c_xs(c), c_ys(c));
        
        subplot(n_rows, n_cols, col);
        plot(b, e_x', linespec{c}, 'color', colors(c, :))
        hold on
        
        subplot(n_rows, n_cols, col+n_cols)
        plot(b, e_y, linespec{c}, 'color', colors(1+c, :)) 
        hold on
    end     
    
    subplot(n_rows, n_cols, col);
    xlabel('$b$ (m)')
    ylabel('RMSE of $\varphi$ (rad)')
    ylim([0 1])
    xticks([-0.5 -0.2 0 0.2 0.5])
    l = legend(c_x_label);
    title(l, '$c_x$')
    subplot(n_rows, n_cols, col+n_cols)
    xlabel('$b$ (m)')
    ylabel('RMSE of $\varphi$ (rad)')
    ylim([0 1])
    xticks([-0.5 -0.2 0 0.2 0.5])
    l = legend(c_y_label);
    title(l, '$c_y$')
end

add_subplot_grid_labels(n_rows, n_cols, {'$v_x$', '$v_y$'}, distance_str)

print('orientation_errors.jpg', '-djpeg','-r600')

function [b, e_x, e_y] = phi_error(values, c_x, c_y)
    
    phi_x = mod(atan2(c_x .* values.x_odd, values.x_even), 2*pi); 
    phi_y = mod(atan2(c_y .* values.y_nav, values.y_odd), 2*pi);
    d_phi_x = compute_angle_difference(values.phi, phi_x);
    d_phi_y = compute_angle_difference(values.phi, phi_y);
    [g, b] = findgroups(values.b);
    e_x = splitapply(@rms, d_phi_x, g);
    e_y = splitapply(@rms, d_phi_y, g);
end
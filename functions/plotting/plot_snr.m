function plot_snr(snr, params, conf)
%plot_snr
%
% Syntax:  plot_snr(snr, params, conf)

%% Easy access variables
sensors = params.sensors;
domain = params.domain;

x_values = unique(snr.sources_bin(:, 1));
y_values = unique(snr.sources_bin(:, 2));
azi_values = unique(snr.sources_bin(:, 3));
cell_size = conf.cell_size;
x_values = x_values(1):cell_size:x_values(end);
y_values = y_values(1):cell_size:y_values(end);
azi_values = azi_values(1):cell_size*pi:azi_values(end);

n_y = length(y_values);
n_x = length(x_values);
n_a = length(azi_values);

%% Compute median in x-y projection
% Determine which cells are empty
[X, Y1] = meshgrid(x_values, y_values);
all_xy_locations = [X(:) Y1(:)];
idx1 = ismembertol(all_xy_locations, snr.sources_bin(:, [1,2]), 0.001, 'ByRows', true);
to_add_xy = all_xy_locations(~idx1,:);
nans_xy = nan(size(to_add_xy, 1), 8);
if any(~idx1)
    warning('not all locations are covered')
end

% Fill the xy table
ratios_xy = table();
ratios_xy.snr_x = snr.x(:, :);
ratios_xy.snr_y = snr.y(:, :);
ratios_xy.x = snr.sources_bin(:,1);
ratios_xy.y = snr.sources_bin(:,2);

ratios_xy = [...
    ratios_xy;...
    table(nans_xy, nans_xy, to_add_xy(:, 1), to_add_xy(:, 2),...
          'VariableNames', {'snr_x', 'snr_y', 'x', 'y'})...
    ];
    
ids = findgroups(ratios_xy.x, ratios_xy.y);
mymedian = @(x) median(x(:), 'omitnan');
values = splitapply(mymedian, ratios_xy.snr_x, ids);
xy_median_x = reshape(values, [n_y, n_x]);
values = splitapply(mymedian, ratios_xy.snr_y, ids);
xy_median_y = reshape(values, [n_y, n_x]);

%% Compute median in azi-y projection
% Determine which cells are empty
[AZI, Y2] = meshgrid(azi_values, y_values);
all_aziy_locations = [AZI(:), Y2(:)];
idx2 = ismembertol(all_aziy_locations, snr.sources_bin(:, [3,2]), 0.001, 'ByRows', true);
to_add_aziy = all_aziy_locations(~idx2,:);
nans_aziy = nan(size(to_add_aziy, 1), 8);
if any(~idx2)
    warning('not all locations are covered')
end

% Filter
% idx = snr.sources_bin(:,1) > -0.2 & snr.sources_bin(:, 1) < 0.2;
idx = 1:size(snr.sources_bin, 1);

% Fill the aziy table
ratios_aziy = table();
ratios_aziy.snr_x = snr.x(idx, :);
ratios_aziy.snr_y = snr.y(idx, :);
ratios_aziy.azi = snr.sources_bin(idx,3);
ratios_aziy.y = snr.sources_bin(idx,2);

ratios_aziy = [...
    ratios_aziy;...
    table(nans_aziy, nans_aziy, to_add_aziy(:, 1), to_add_aziy(:, 2),...
          'VariableNames', {'snr_x', 'snr_y', 'azi', 'y'})...
    ];
    
ids = findgroups(ratios_aziy.azi, ratios_aziy.y);
mymedian = @(x) median(x(:), 'omitnan');
values = splitapply(mymedian, ratios_aziy.snr_x, ids);
aziy_median_x = reshape(values, [n_y, n_a]);
values = splitapply(mymedian, ratios_aziy.snr_y, ids);
aziy_median_y = reshape(values, [n_y, n_a]);

%% Plot x-y projection

% Correct for sensor location
y_range = domain.y_range - sensors.y_value;

% Make x figure
figure
imagesc(x_values, y_values - sensors.y_value, xy_median_x);
hold on
% contour(X, Y1, xy_median_x, [1 1], 'w-')
colormap(conf.map)
set(gca, 'ydir', 'normal')
axis image
hold on
mk = plot(sensors.locations.x, 0, 'ko', 'markerfacecolor', [0 0 0]);
xlim([min(domain.x_range), max(domain.x_range)])
ylim([-0.025 y_range(2)+0.025])
xlabel('$x$~(m)') 
ylabel('$y$~(m)')
caxis(conf.color_range)   
set(gca, 'tickdir', 'out')
    
drawnow
post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
curunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
cursize = get(gca, 'Position');
set(gca, 'Units', curunits);
pt_sz = diff(xlim) / cursize(3);
sensor_sz = 8e-3 / pt_sz;
set(mk, 'MarkerSize', sensor_sz)
% xticks(-2:0.2:2)
% yticks(0:0.2:1)
name = ['snr_spatial_x.jpg'];
print(fullfile(params.output_folder, name), '-djpeg', '-r600')
disp(name)
close

% Make y figure
figure
imagesc(x_values, y_values - sensors.y_value, xy_median_y);
hold on
% contour(X, Y1, xy_median_x, [1 1], 'w-')
colormap(conf.map)
set(gca, 'ydir', 'normal')
axis image
hold on
mk = plot(sensors.locations.x, 0, 'ko', 'markerfacecolor', [0 0 0]);
xlim([min(domain.x_range), max(domain.x_range)])
ylim([-0.025 y_range(2)+0.025])
xlabel('$x$~(m)') 
ylabel('$y$~(m)')
caxis(conf.color_range)
set(gca, 'tickdir', 'out')

drawnow
post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
curunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
cursize = get(gca, 'Position');
set(gca, 'Units', curunits);
pt_sz = diff(xlim) / cursize(3);
sensor_sz = 8e-3 / pt_sz;
set(mk, 'MarkerSize', sensor_sz)
% xticks(-2:0.2:2);
% yticks(0:0.2:1)
name = ['snr_spatial_y.jpg'];
print(fullfile(params.output_folder, name), '-djpeg', '-r600')
disp(name)
% close

% Make colorbar
figure
colormap(conf.map);
caxis(conf.color_range)
c = colorbar('north');
c.Label.Interpreter = 'Latex';
c.Label.String = 'median SNR~(dB)';
axis off
drawnow
post_process_figure(0.8, 0.11, [0 0], [0 -1.5])
name = ['snr_colorbar.jpg'];
print(fullfile(params.output_folder, name), '-djpeg', '-r600')
disp(name)
close

%% Plot azi-y projection

% Convert azi values to degrees and shift it for proper orientation
azi_in_degrees = -azi_values / pi * 180 +  90;

% Make x figure
figure
polarPcolor(y_values-sensors.y_value, azi_in_degrees, aziy_median_x, ...
    'colormap', conf.map,...
    'NSpokes', conf.n_spokes,...
    'NCircles', conf.n_circles,...
    'colBar', 0);
caxis(conf.color_range)
post_process_figure(0.4, 0.65, [0 0], [0.4 -0.3])
text(-0.15, 1.26, '$\varphi$ (rad)')
text(0.3, -0.1', '$y - s_y$ (m)')
name = ['snr_orientation_x.jpg'];
print(fullfile(params.output_folder, name), '-djpeg', '-r600')
disp(name)
close

% Make y figure
figure
polarPcolor(y_values-sensors.y_value, azi_in_degrees, aziy_median_y, ...
    'colormap', conf.map,...
    'NSpokes', conf.n_spokes,...
    'NCircles', conf.n_circles,...
    'colBar', 0);
caxis(conf.color_range)
post_process_figure(0.4, 0.65, [0 0], [0.4 -0.3])
text(-0.15, 1.26, '$\varphi$ (rad)')
text(0.3, -0.1', '$y - s_y$ (m)')
name = ['snr_orientation_y.jpg'];
print(fullfile(params.output_folder, name), '-djpeg', '-r600')
disp(name)
close

end
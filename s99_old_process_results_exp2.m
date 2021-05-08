%% Introduction
% Create plots of the input mode experiment

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

% Input files
folder = './data/20_05_21/';

% Params
load([folder, 'params_input_mode_1.mat'], 'params');

% Plotting variables
cell_size = 0.02;
ranges = [1 3 5 9]; % in source, 1 source outside, 2 source outside, 4 source outside
median_error_thresholds = [0 params.source.radius .* ranges 1];
error_score_thresholds = [0 0.05 0.25 0.5 0.75 1];
map = magma(1024);


%% Merge output test_methods_data
files = dir([folder, 'predictions_input_mode*.mat']);
merged = [];
for idx = 1:length(files)
    load([folder, files(idx).name])
    name = files(idx).name;
    mode = round(str2double(name(end-4:end-4)));
    predictions.input_mode = repmat(params.experiment.input_modes(mode), height(predictions), 1);
    merged = [merged; predictions];
    clear predictions
end
predictions = merged;
predictions.Properties.VariableNames{1} = 'method';
predictions.method = categorical(predictions.method);
predictions.input_mode = categorical(predictions.input_mode);

% Add error and difference
predictions.error = compute_prediction_errors(predictions.prediction, predictions.source, params);
predictions.difference = predictions.prediction - predictions.source;
predictions.difference(:,3) = compute_angle_difference(predictions.source(:, 3), predictions.prediction(:, 3)) / pi;
predictions.component_error = abs(predictions.difference);

% Add location bins
predictions.source_bin = compute_bins(predictions.source, cell_size);

[~, methods, input_mode] = findgroups(predictions.method, predictions.input_mode);
id_cat = arrayfun(@(x) [char(methods(x)) ' ' char(input_mode(x))], 1:length(methods), 'UniformOutput', false);
id_cat = categorical(id_cat);

%% Median error plots
% Compute spatial median errors
% TODO requires that all combinations are present in predictions!
[xy_bins, method, input_mode, x, y] = findgroups(...
    predictions.method,...
    predictions.input_mode,...
    predictions.source_bin(:,1),...
    predictions.source_bin(:,2));
median_error = splitapply(@median, predictions.error, xy_bins);
spatial_error = table(method, input_mode, x, y, median_error);

filter = predictions.source(:,1) > -0.2 & predictions.source(:,1) < 0.2;
[ay_bins, method, input_mode, azi, y] = findgroups(...
    predictions.method(filter,:),...
    predictions.input_mode(filter,:),...
    predictions.source_bin(filter,3),...
    predictions.source_bin(filter,2));
median_error = splitapply(@median, predictions.error(filter,:), ay_bins);
orientational_error = table(method, input_mode, azi, y, median_error);

x = unique(x);
y = unique(y);
a = unique(azi);

% % spatial heatmap plot
% heatmap_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_spatial_heatmap_input_mode_' num2str(mode) '.jpg'];
% conf.caxis = [0 median_error_thresholds(end-1)];
% conf.colormap = flip(map);
% conf.cell_size = cell_size;
% 
% [ids, methods, input_modes] = findgroups(spatial_error.method, spatial_error.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = heatmap_name(char(method), mode_idx);
%    plot_spatial_heatmap(x, y, spatial_error.median_error(filter), params, conf);
% end
% 
% % orientation heatmap plot
% orientation_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_orientation_heatmap_input_mode_' num2str(mode) '.jpg'];
% conf.caxis = [0 median_error_thresholds(end - 1)];
% conf.colormaps = flip(map);
% conf.cell_size = conf.cell_size;
% conf.spokes = 5;
% conf.circles = 3;
% 
% [ids, methods, input_modes] = findgroups(orientational_error.method, orientational_error.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = orientation_name(char(method), mode_idx);
%    plot_orientation_heatmap(a, y, orientational_error.median_error(filter), params, conf);
% end
% 
% % heatmap colorbar
% figure
% colormap(flip(map));
% caxis([0 median_error_thresholds(end-1)])
% c = colorbar('north');
% c.Label.Interpreter = 'Latex';
% c.Label.String = 'median error';
% c.Label.FontSize = 9;
% axis off
% drawnow
% post_process_figure(0.8, 0.12, [0 0], [0 -1.5])
% drawnow
% pause(0.2)
% print('./images/exp2_input_mode/error_heatmap_colorbar.jpg', '-djpeg', '-r600')
% disp('error_heatmap_colorbar')
% close

% spatial contour plot
contour_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_spatial_contour_input_mode_' num2str(mode) '.jpg'];
conf.thresholds = median_error_thresholds;
conf.cell_size = cell_size;
conf.colormap = flip(map);

[ids, methods, input_modes] = findgroups(spatial_error.method, spatial_error.input_mode);
u_id = unique(ids);
for id = u_id'
   filter = ids == id; 
   method = methods(id);
   input_mode = input_modes(id);
   mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
   conf.name = contour_name(char(method), mode_idx);
   plot_spatial_contour(x, y, spatial_error.median_error(filter), params, conf);
end

% orientation contour plot
orientation_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_orientation_contour_input_mode_' num2str(mode) '.jpg'];
conf.caxis = [0 median_error_thresholds(end - 1)];
conf.colormaps = flip(map);
conf.cell_size = conf.cell_size;
conf.spokes = 5;
conf.circles = 3;

[ids, methods, input_modes] = findgroups(orientational_error.method, orientational_error.input_mode);
u_id = unique(ids);
for id = u_id'
   filter = ids == id; 
   method = methods(id);
   input_mode = input_modes(id);
   mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
   conf.name = orientation_name(char(method), mode_idx);
   plot_orientation_contour(a, y, orientational_error.median_error(filter), params, conf);
end

% contour colorbar
my_map = flip(map);
cmax = median_error_thresholds(end-1);
cmin = 0;
index = fix((median_error_thresholds - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
colors = squeeze(ind2rgb(index', my_map));
colors = [1 1 1; colors];
d_colors = get(groot, 'defaultAxesColorOrder');

figure
values = median_error_thresholds - [0 median_error_thresholds(1:end-1)];
set(groot, 'defaultAxesColorOrder', colors); 
barh(repmat(values, 2, 1), 0.05, 'stacked')
set(groot, 'defaultAxesColorOrder', d_colors);
set(gca, 'tickdir', 'in')
axis image
ylim([0.98, 1.025]);
xlim([0, median_error_thresholds(end-1)+0.11]);
yticks([])
xticks([median_error_thresholds(1:end-1) max(xlim)])
xlabel('median error')
daspect([1 8.5 1])
post_process_figure(0.8, 0.12, [0.3 0.65], [0.3 0])

drawnow
pause(0.2)
print('./images/exp2_input_mode/error_contour_colorbar.jpg', '-djpeg', '-r600')
disp('error_contour_colorbar')
close

figure
values = median_error_thresholds - [0 median_error_thresholds(1:end-1)];
values(end) = [];
set(groot, 'defaultAxesColorOrder', colors); 
barh(repmat(values, 2, 1), 0.05, 'stacked')
set(groot, 'defaultAxesColorOrder', d_colors);
set(gca, 'tickdir', 'in')
axis image
ylim([0.98, 1.025]);
xlim([0, median_error_thresholds(end-1)]);
yticks([])
xticks(median_error_thresholds(1:end-1))
xlabel('median error')
daspect([1 8.5 1])
post_process_figure(0.4, 0.24, [0.3 0.65], [0.3 0])

drawnow
pause(0.2)
print('./images/exp2_input_mode/error_bar_colorbar.jpg', '-djpeg', '-r600')
disp('error_bar_colorbar')
close
% % spatial block plot
% block_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_spatial_block_input_mode_' num2str(mode) '.jpg'];
% conf.thresholds = median_error_thresholds;
% conf.colormap = map;
% conf.cell_size = cell_size; 
% 
% [ids, methods, input_modes] = findgroups(spatial_error.method, spatial_error.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = block_name(char(method), mode_idx);
%    plot_spatial_block(x, y, spatial_error.median_error(filter), params, conf);
% end
% 
% % orientation block plot
% block_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_error_orientation_block_input_mode_' num2str(mode) '.jpg'];
% conf.thresholds = median_error_thresholds;
% conf.colormap = map;
% conf.cell_size = cell_size; 
% 
% [ids, methods, input_modes] = findgroups(orientational_error.method, orientational_error.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = block_name(char(method), mode_idx);
%    plot_orientation_block(unique(a), unique(y), orientational_error.median_error(filter), params, conf);
% end

% Area error-thresholds
thresh = [median_error_thresholds(1:end-1) Inf];
[ids, method, input_mode] = findgroups(spatial_error.method, spatial_error.input_mode);
area = splitapply(@(x) sum(x <= thresh), spatial_error.median_error, ids);
n_groups = length(unique(ids));
num_cells = height(spatial_error) / n_groups;
max_area = num_cells * cell_size^2;
area = area - [zeros(n_groups, 1) area(:, 1:end-1)];
area = area * cell_size^2; % ./ max_area * 100;
method(method == 'rand') = 'rnd';
threshold_areas = table(method, input_mode, area);

% Area bar plot
my_map = flip(map);
cmax = median_error_thresholds(end-1);
cmin = 0;
index = fix((median_error_thresholds - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
colors = squeeze(ind2rgb(index', my_map));
colors = [1 1 1; colors];
d_colors = get(groot, 'defaultAxesColorOrder');

threshold_areas.method = upper(arrayfun(@char, threshold_areas.method, 'UniformOutput', false));
threshold_areas.method = categorical(threshold_areas.method);
u_methods = unique(threshold_areas.method);
n_methods = length(u_methods);
u_modes = unique(threshold_areas.input_mode);
n_modes = length(u_modes);

method_order = {'GN', 'LSQ', 'NR', 'KNN', 'MLP', 'ELM', 'LCMV', 'CWT', 'RND'};
assert(length(method_order) == n_methods);
method_order = arrayfun(@(idx)...
    find(strcmp(method_order{idx},...
                arrayfun(@char, u_methods, 'UniformOutput', false))),...
    1:length(method_order));
mode_characters = {'x', 'b', 'a', 'y'};
mode_order = [2 1 4 3];

area = permute(threshold_areas.area, [1 3 2]);
area = reshape(area, [n_modes, n_methods, length(median_error_thresholds)]);
area = area(mode_order, method_order, :);
area = permute(area, [2 1 3]);

area(:, :, end) = [];

method_labels = arrayfun(@(x) {char(x)}, u_methods(method_order));

mode_labels = mode_characters(mode_order);
%arrayfun(@(x) {num2str(find(strcmp(params.experiment.input_modes, char(x))))}, u_modes(mode_order,:));

set(groot, 'defaultAxesColorOrder', colors); 
figure
plotBarStackGroups(area, method_labels, mode_labels);
ylim([0 max_area])
ylabel('area (m$^2$)')
xlabel('method / input mode')
post_process_figure(1, 0.45, [1.2 1.4], [0 0.2])
set(groot, 'defaultAxesColorOrder', d_colors);

drawnow
pause(0.2)
print('./images/exp2_input_mode/median_error_area_bars', '-djpeg', '-r600')
disp('median_error_area_bars')
close

% %% Error score plots
% % Compute spatial scores
% % TODO requires that all combinations are present in predictions!
% [xy_bins, method, input_mode, x, y] = findgroups(...
%     predictions.method,...
%     predictions.input_mode,...
%     predictions.source_bin(:,1),...
%     predictions.source_bin(:,2));
% score = splitapply(@(x) sum(x > params.source.radius) / length(x), predictions.error, xy_bins);
% spatial_score = table(method, input_mode, x, y, score);
% 
% % Compute orientation scores
% filter = predictions.source(:,1) > -0.2 & predictions.source(:,1) < 0.2;
% [ay_bins, method, input_mode, azi, y] = findgroups(...
%     predictions.method(filter,:),...
%     predictions.input_mode(filter,:),...
%     predictions.source_bin(filter,3),...
%     predictions.source_bin(filter,2));
% score = splitapply(@(x) sum(x > params.source.radius) / length(x), predictions.error(filter,:), ay_bins);
% orientational_score = table(method, input_mode, azi, y, score);
% 
% x = unique(x);
% y = unique(y);
% a = unique(azi);
% 
% % Spatial heatmap
% heatmap_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_score_spatial_heatmap_input_mode_' num2str(mode) '.jpg'];
% conf.caxis = [0 1];
% conf.colormap = flip(map);
% conf.cell_size = cell_size;
% 
% [ids, methods, input_modes] = findgroups(spatial_score.method, spatial_score.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = heatmap_name(char(method), mode_idx);
%    plot_spatial_heatmap(x, y, spatial_score.score(filter), params, conf);
% end
% 
% % orientation  heatmap
% orientation_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_score_orientation_heatmap_input_mode_' num2str(mode) '.jpg'];
% conf.caxis = [0 1];
% conf.colormaps = flip(map);
% conf.cell_size = conf.cell_size;
% conf.spokes = 5;
% conf.circles = 3;
% 
% [ids, methods, input_modes] = findgroups(orientational_score.method, orientational_score.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = orientation_name(char(method), mode_idx);
%    plot_orientation_heatmap(a, y, orientational_score.score(filter), params, conf);
% end
% 
% % heatmap colorbar
% figure
% colormap(flip(map));
% caxis([0 1])
% c = colorbar('north');
% c.Label.Interpreter = 'Latex';
% c.Label.String = 'error score';
% c.Label.FontSize = 9;
% axis off
% drawnow
% post_process_figure(0.8, 0.12, [0 0], [0 -1.5])
% drawnow
% pause(0.2)
% print('./images/exp2_input_mode/score_heatmap_colorbar.jpg', '-djpeg', '-r600')
% disp('score_heatmap_colorbar')
% close
% 
% % spatial contour
% contour_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_score_spatial_contour_input_mode_' num2str(mode) '.jpg'];
% conf.thresholds = error_score_thresholds;
% conf.cell_size = cell_size;
% conf.colormap = flip(map);
% 
% [ids, methods, input_modes] = findgroups(spatial_score.method, spatial_score.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = contour_name(char(method), mode_idx);
%    plot_spatial_contour(x, y, spatial_score.score(filter), params, conf);
% end
% 
% % orientation contour
% orientation_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_score_orientation_contour_input_mode_' num2str(mode) '.jpg'];
% conf.caxis = [0 error_score_thresholds(end - 1)];
% conf.colormaps = flip(map);
% conf.cell_size = cell_size;
% conf.spokes = 5;
% conf.circles = 3;
% 
% [ids, methods, input_modes] = findgroups(orientational_score.method, orientational_score.input_mode);
% u_id = unique(ids);
% for id = u_id'
%    filter = ids == id; 
%    method = methods(id);
%    input_mode = input_modes(id);
%    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
%    conf.name = orientation_name(char(method), mode_idx);
%    plot_orientation_contour(a, y, orientational_score.score(filter), params, conf);
% end
% 
% % contour colorbar
% my_map = flip(map);
% cmax = error_score_thresholds(end-1);
% cmin = 0;
% index = fix((error_score_thresholds - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
% colors = squeeze(ind2rgb(index', my_map));
% d_colors = get(groot, 'defaultAxesColorOrder');
% set(groot, 'defaultAxesColorOrder', colors); 
% 
% figure
% values = error_score_thresholds - [0 error_score_thresholds(1:end-1)];
% barh(repmat(values(2:end), 2, 1), 0.05, 'stacked')
% set(gca, 'tickdir', 'in')
% axis image
% ylim([0.98, 1.025]);
% xlim([0, 1]);
% yticks([])
% xticks([0 error_score_thresholds(2:end-1) 1])
% xlabel('error score')
% daspect([1 1.65 1])
% post_process_figure(0.8, 0.12, [0.3 0.65], [0.3 -0.15])
% set(groot, 'defaultAxesColorOrder', d_colors);
% 
% drawnow
% pause(0.2)
% print('./images/exp2_input_mode/score_contour_colorbar.jpg', '-djpeg', '-r600')
% disp('error_contour_colorbar')
% close
% 
% % Area score-thresholds
% thresh = error_score_thresholds;
% [ids, method, input_mode] = findgroups(spatial_score.method, spatial_score.input_mode);
% area = splitapply(@(x) sum(x <= thresh), spatial_score.score, ids);
% n_groups = length(unique(ids));
% num_cells = height(spatial_score) / n_groups;
% max_area = num_cells * cell_size^2;
% area = area - [zeros(n_groups, 1) area(:, 1:end-1)];
% area = area * cell_size^2 ./ max_area * 100;
% threshold_areas = table(method, input_mode, area);
% 
% % Area score bar plot
% my_map = flip(map);
% cmax = error_score_thresholds(end-1);
% cmin = 0;
% index = fix((error_score_thresholds(1:end) - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
% colors = squeeze(ind2rgb(index', my_map));
% d_colors = get(groot, 'defaultAxesColorOrder');
% 
% u_methods = unique(threshold_areas.method);
% n_methods = length(u_methods);
% u_modes = unique(threshold_areas.input_mode);
% n_modes = length(u_modes);
% 
% method_order = 1:n_methods; %[4 3 6 5 7 2 1 8];
% mode_order = [2 1 4 3];
% 
% area = permute(threshold_areas.area, [1 3 2]);
% area(:, :, 1) = area(:, :, 2) + area(:, :, 1);
% area(:, :, 2) = [];
% area = reshape(area, [n_modes, n_methods, length(error_score_thresholds)-1]);
% area = area(mode_order, method_order, :);
% area = permute(area, [2 1 3]);
% 
% method_labels = arrayfun(@(x) {char(x)}, u_methods(method_order));
% mode_labels = arrayfun(@(x) {num2str(find(strcmp(params.experiment.input_modes, char(x))))}, u_modes(mode_order,:));
% 
% set(groot, 'defaultAxesColorOrder', colors); 
% figure
% plotBarStackGroups(area, method_labels, mode_labels);
% ylim([0 35])
% ylabel('area (\% of total)')
% xlabel('method / input mode')
% post_process_figure(1, 0.5, [1 1.4], [0 0.2])
% set(groot, 'defaultAxesColorOrder', d_colors);
% 
% drawnow
% pause(0.2)
% print('./images/exp2_input_mode/error_score_area_bars', '-djpeg', '-r600')
% disp('error_score_area_bars')
% close

%% Predicted locations
num_locations = 200;
prediction_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_predictions_input_mode_' num2str(mode) '.jpg'];
[ids, methods, input_modes] = findgroups(predictions.method, predictions.input_mode);
u_id = unique(ids);
max_index = sum(ids == 1);
show_filter = randsample(max_index, num_locations);
for id = u_id'
    filter = find(ids == id);
    assert(length(filter) == max_index);
    filter = filter(show_filter);
    input_mode = input_modes(id);
    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
    method = methods(id);
    plot_predicted_locations(predictions.prediction(filter,:), predictions.source(filter,:), params)
    name = prediction_name(char(method), mode_idx);
    
    drawnow
    pause(0.2)
    print(name, '-djpeg', '-r600')
    disp(name)
    close
end


%% Plot location and orientation error individually
% Spatial medians
location_error = sqrt(sum((predictions.prediction(:, 1:2) - predictions.source(:, 1:2)).^2, 2));
orientation_error = abs(compute_angle_difference(predictions.source(:, 3), predictions.prediction(:, 3)) / pi);

[xy_bins, method, input_mode, x, y] = findgroups(...
    predictions.method,...
    predictions.input_mode,...
    predictions.source_bin(:,1),...
    predictions.source_bin(:,2));
median_location_error = splitapply(@median, location_error, xy_bins);
median_orientation_error = splitapply(@median, orientation_error, xy_bins);
spatial_errors = table(method, input_mode, x, y, median_location_error, median_orientation_error);

% Define colormap
conf.location_error_min = 0;
conf.location_error_max = median_error_thresholds(end - 1);
conf.orientation_error_min = 0;
conf.orientation_error_max = median_error_thresholds(end - 1);
conf.color_low_low = [1 1 1];    % RGB
conf.color_high_high = [0 0 0];  % RGB
conf.color_low_high = [0.8 0 0];   % RGB
conf.color_high_low = [0 0 0.8];   % RGB
conf.cell_size = cell_size;
conf.thresholds = median_error_thresholds;

% Compute threshold combinations
[L, O] = meshgrid(median_error_thresholds, median_error_thresholds);
% Define base colors for interpolation
[X, Y] = meshgrid([conf.location_error_min conf.location_error_max], [conf.orientation_error_min conf.orientation_error_max]);
v = [conf.color_low_low; conf.color_low_high; conf.color_high_low; conf.color_high_high];
% Interpolate other values
color_comb = [interp2(X, Y, reshape(v(:, 1), [2 2]) , L(:), O(:), 'linear', v(end, 1)),...
              interp2(X, Y, reshape(v(:, 2), [2 2]) , L(:), O(:), 'linear', v(end, 2)),...
              interp2(X, Y, reshape(v(:, 3), [2 2]) , L(:), O(:), 'linear', v(end, 3))];
n_colors = length(median_error_thresholds);

% Plot colormap
pcolor(L, O, reshape(1:size(color_comb, 1), [n_colors, n_colors]))
colormap(color_comb)
set(gca, 'ydir', 'normal')
xlabel('location error (m)')
ylabel('orientation error ($\times\pi$ rad)')
xlim([0 0.2])
ylim([0 0.2])
ticks = median_error_thresholds;
ticks(end) = 0.2; ticks(1) = [];
xticks(ticks)
yticks(ticks)
post_process_figure(0.24, 1, [0.85 0.6], [0.2 0])
set(gca, 'tickdir', 'out')
drawnow
pause(0.2)
name = './images/exp2_input_mode/colorcoded_colorbar.jpg';
print(fullfile(name), '-djpeg', '-r600')
disp(name)
close

% Helper values
u_x = unique(x);
u_y = unique(y);
n_x = length(u_x);
n_y = length(u_y);

% Make the plots for each method and resolution
colorcoded_name = @(mtd, mode) ['./images/exp2_input_mode/' lower(mtd) '_colorcoded_spatial_input_mode_' num2str(mode) '.jpg'];
[ids, methods, input_modes] = findgroups(spatial_errors.method, spatial_errors.input_mode);
u_id = unique(ids);
for id = u_id'
    filter = ids == id; 
    method = methods(id);
    input_mode = input_modes(id);
    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
    conf.name = colorcoded_name(char(method), mode_idx);
   
    L = spatial_errors.median_location_error(filter);
    O = spatial_errors.median_orientation_error(filter);
   
    color_index = zeros(size(L, 1), 1);
%     loc_lvls = zeros(size(L, 1), 1);
%     or_lvls = zeros(size(L, 1), 1);
    for idx = 1:size(L, 1)
       loc_lvl = find(L(idx, :) <= median_error_thresholds, 1, 'first') - 1;
       or_lvl =  find(O(idx, :) <= median_error_thresholds, 1, 'first') - 1;
       if isempty(loc_lvl)
           loc_lvl = 6;
       end
       if isempty(or_lvl)
           or_lvl = 6;
       end
       
       color_index(idx) = sub2ind([n_colors, n_colors], or_lvl, loc_lvl);
%        loc_lvls(idx) = loc_lvl;
%        or_lvls(idx) = or_lvl;
    end
    color_index = reshape(color_index, [n_y, n_x]);
    figure
    imagesc(u_x, u_y, color_index)
    colormap(color_comb)
    caxis([1 size(color_comb, 1)])
    set(gca, 'ydir', 'normal')
    hold on
    plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
    ylim([0-conf.cell_size / 2 max(y) + conf.cell_size / 2]);
    xlabel('$x$ (m)')
    ylabel('$y$ (m)')
    daspect([1 1 1])
    set(gca, 'tickdir', 'out')
    post_process_figure(0.24, 0.63, [0.8 0.6], [0.15 0])
    drawnow
    pause(0.2)
    print(fullfile(conf.name), '-djpeg', '-r600')
    disp(conf.name)
    close
end

% Orientation median
filter = predictions.source(:,1) > -0.2 & predictions.source(:,1) < 0.2;
[ay_bins, method, input_mode, azi, y] = findgroups(...
    predictions.method(filter,:),...
    predictions.input_mode(filter,:),...
    predictions.source_bin(filter,3),...
    predictions.source_bin(filter,2));
median_location_error = splitapply(@median, location_error(filter,:), ay_bins);
median_orientation_error = splitapply(@median, orientation_error(filter, :), ay_bins);
orientation_errors = table(method, input_mode, azi, y, median_location_error, median_orientation_error);

% Helper values
u_a = unique(azi);
u_y = unique(y);
n_a = length(u_a);
n_y = length(u_y);

conf.colormap = color_comb;
conf.cell_size = cell_size;
conf.spokes = 5;
conf.circles = 3;

% Make the plots for each method and resolution
colorcoded_name = @(mtd, idx) ['./images/exp2_input_mode/' lower(mtd) '_colorcoded_orientation_input_mode_' num2str(idx) '.jpg'];
[ids, methods, input_modes] = findgroups(orientation_errors.method, orientation_errors.input_mode);
u_id = unique(ids);
for id = u_id'
    filter = ids == id; 
    method = methods(id);
    input_mode = input_modes(id);
    mode_idx = find(strcmp(params.experiment.input_modes, char(input_mode)));
    conf.name = colorcoded_name(char(method), mode_idx);
   
    L = orientation_errors.median_location_error(filter);
    O = orientation_errors.median_orientation_error(filter);
    
    color_index = zeros(size(L, 1), 1);
    for idx = 1:size(L, 1)
       loc_lvl = find(L(idx, :) <= median_error_thresholds, 1, 'first') - 1;
       or_lvl =  find(O(idx, :) <= median_error_thresholds, 1, 'first') - 1;
       if isempty(loc_lvl)
           loc_lvl = 6;
       end
       if isempty(or_lvl)
           or_lvl = 6;
       end
       
       color_index(idx) = sub2ind([n_colors, n_colors], or_lvl, loc_lvl);
    end
    color_index = reshape(color_index, [n_y, n_a]);
    
    plot_orientation_colorcoded(u_a, u_y, color_index, params, conf);
end

%% Compute error distributions
% values = splitapply(@(x) prctile(x, [5 25 50 75 95]), predictions.error, ids);
% 
% distributions = table(...
%     methods, input_mode, values(:, 1), values(:, 2),...
%     values(:, 3), values(:, 4), values(:, 5));
% distributions.Properties.VariableNames = {'method' 'resolution' 'p_5' 'p_25' 'p_50' 'p_75' 'p_95'};
% 
% figure
% bar(id_cat, distributions{:, 3:end}, 'stacked')

%% 

reduced = predictions(:, {'method', 'input_mode', 'error'});
writetable(reduced, 'reduced.csv')
function handles = violin(data, conf)
%violin
%
% Syntax: handles = violin(data, group_labels, violin_labels, conf)
%
% Creates grouped violin plots
%
% - Updated by Jelmer Bot (2019)
% Adapted from:
%__________________________________________________________________________
% violin.m - Simple violin plot using matlab default kernel density estimation
% Last update: 10/2015
%__________________________________________________________________________
% This function creates violin plots based on kernel density estimation
% using ksdensity with default settings. Please be careful when comparing pdfs
% estimated with different bandwidth!
%__________________________________________________________________________
%
% Please cite this function as:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%__________________________________________________________________________
%
%% Configurable stuff
violin_width = 0.8; % 1 is touching
group_width = 0.85;  % 1 is touching

%% Internal values
num_groups = size(data, 1);
num_violins_per_group = size(data, 2);
num_violins = num_groups * num_violins_per_group;
conf = process_variables(num_violins, conf);
if length(conf.violin_labels) ~= num_violins_per_group
    error('Then number of violin_labels does not match the number of violins per group');
end
if length(conf.group_labels) ~= num_groups
   error('The number of group_labels does not match the number of groups'); 
end
    
%% Compute density and other values
loc = linspace(conf.y_range(1), conf.y_range(2), conf.num_points)';
smooth_points = round(conf.bandwidth ./ (diff(conf.y_range) ./ conf.num_points));
densities = nan(num_groups, num_violins_per_group, conf.num_points);
% locations = nan(num_groups, num_violins_per_group, conf.num_points);
% bandwidths = nan(num_groups, num_violins_per_group, 1);
% means = nan(num_groups, num_violins_per_group, 1);
% medians = nan(num_groups, num_violins_per_group, 1);
percentiles = nan(num_groups, num_violins_per_group, 5);
for group_idx = 1:num_groups
    for violin_idx = 1:num_violins_per_group
        violin_data = data{group_idx, violin_idx};
%         figure
%         hist(violin_data);
        if ~isempty(violin_data)
%             [density, location, bandwidth] = ...
%                         ksdensity(violin_data,...
%                             'NumPoints', conf.num_points,...
% ...                            'Support', conf.support_range,...
%                             'Bandwidth', conf.bandwidth);
            counts = smooth(histc(violin_data, loc), smooth_points);
            density= counts ./ max(counts);
            densities(group_idx, violin_idx, :) = density;
%             locations(group_idx, violin_idx, :) = loc;
%             bandwidths(group_idx, violin_idx) = bandwidth;
%             means(group_idx, violin_idx) = nanmean(violin_data);
%             medians(group_idx, violin_idx) = nanmedian(violin_data);
            percentiles(group_idx, violin_idx,:) = prctile(violin_data, [5 25 50 75 95]);
        end
    end
end

%% Plot the violins
figure
hold on

violin_positions = zeros(num_violins_per_group, num_groups);

violin_idx = 1;
for group_idx = 1:num_groups
    for violin_in_group_idx = 1:num_violins_per_group
        density = squeeze(densities(group_idx, violin_in_group_idx, :));        
%         mean_value = means(group_idx, violin_in_group_idx);
%         median_value = medians(group_idx, violin_in_group_idx);
        percentile_value = squeeze(percentiles(group_idx, violin_in_group_idx, :)); 
        
        % compute the x location
        total_violin_width = group_width / num_violins_per_group;
        position_in_group = violin_in_group_idx - (num_violins_per_group + 1) / 2;
        x_offset = position_in_group * total_violin_width + group_idx;
        violin_positions(violin_in_group_idx, group_idx) = x_offset;
        
        % Scale the density, skip if information is incomplete
        x_values = density * total_violin_width * violin_width / 2;
        y_values = squeeze(loc);
        if any(isnan(x_values))
            continue
        end
        
        i_start = find(x_values > 0, 1, 'first');
        i_end = find(x_values > 0, 1, 'last');
        x = x_values(i_start:i_end);
        y = y_values(i_start:i_end);
        
        % Plot the violin
        handles.violins(violin_idx) = fill(...
            [-x + x_offset; flipud(x + x_offset)],...
            [y; flipud(y)], ...
            conf.face_color(violin_idx,:),...
            'FaceAlpha', conf.face_alpha);
        
%         density_at_mean = interp1(y_values, x_values, mean_value);
%         handles.means(violin_idx) = plot([...
%             -density_at_mean + x_offset,...
%              density_at_mean+ x_offset...
%             ], repmat(mean_value, 1, 2), conf.mean_colour);
        
        % Plot percentiles
        density_at_percentiles = interp1(y_values, x_values, percentile_value);
        handles.medians(violin_idx, :) = plot([...
            -density_at_percentiles + x_offset,...
             density_at_percentiles + x_offset,...
            ]', repmat(percentile_value, 1, 2)', conf.median_colour);
        handles.medians(violin_idx, 3).LineWidth = 1.2;
        
%         % Plot violin around percentiles
%         x_values = density_at_percentiles;
%         y_values = percentile_value;
%         handles.violines(violin_idx) = fill(...
%             [-x_values + x_offset; flipud(x_values + x_offset)],...
%             [y_values; flipud(y_values)], ...
%             conf.face_color(violin_idx,:),...
%             'FaceAlpha', conf.face_alpha);
%         
        violin_idx = violin_idx + 1;
        
    end
end

%% Create the labels and ticks
half_group_size = group_width / 2;
spacing = 1-group_width;
xlim([1 - half_group_size - spacing, num_groups + half_group_size + spacing]); 
ax = gca;
ax.TickDir = 'out';                % Stylistic choice
ax.TickLabelInterpreter = 'Latex'; % Required for multi-line tick-labels

make_label = @(l1, l2) ['\begin{tabular}{c} ' l1 ' \\ ' l2 ' \end{tabular}'];

if mod(num_violins_per_group, 2) == 1 % isodd 
    % The center tick, which gets the group-label is in
    % violin_positions. So we only need to prepare the labels
    half = floor(num_violins_per_group / 2);
    center = half + 1;
    labels = cell(1, num_groups * num_violins_per_group);
    
    for idx = 1:length(labels)
        violin_idx = mod(idx - 1, num_violins_per_group) + 1;
        group_idx = floor((idx - 1) / num_violins_per_group) + 1;
        
        if violin_idx == center
            labels(idx) = {make_label(conf.violin_labels{violin_idx}, conf.group_labels{group_idx})};
        else
            labels(idx) = {make_label(conf.violin_labels{violin_idx}, ' ')};
        end
    end
    
    ax.XTick = violin_positions(:);
    ax.XTickLabels = labels;
else
    % The group-tick is placed between violins, so add it manually
    num_ticks = num_groups * (num_violins_per_group + 1);
    positions = zeros(1, num_ticks);
    half = floor(num_violins_per_group / 2);
    center = half + 1;
    group_tick_idx = center:(num_violins_per_group + 1):num_ticks;
    positions(group_tick_idx) = 1:num_groups;
    positions(positions == 0) = violin_positions(:);
    
    % Prepare the labels
    labels = cell(1, num_groups * (num_violins_per_group + 1));
    for idx = 1:length(labels)
        violin_idx = mod(idx - 1, (num_violins_per_group + 1)) + 1;
        group_idx = floor((idx - 1) / (num_violins_per_group + 1)) + 1;
        
        if violin_idx == center
            % The center tick for the group label
            labels(idx) = {make_label(' ', conf.group_labels{group_idx})};
        elseif violin_idx < center
            % Base case
            labels(idx) = {make_label(conf.violin_labels{violin_idx}, ' ')};
        else
            % We insterted a tick manually, so correct the violin_index
            labels(idx) = {make_label(conf.violin_labels{violin_idx - 1}, ' ')};
        end
    end
    
    ax.XTick = positions;
    ax.XTickLabels = labels;
end

end

function conf = process_variables(num_violins, conf)

% Default values
if ~isfield(conf, 'face_color')
    conf.face_color=[1 0.5 0];
end
if ~isfield(conf, 'face_alpha')
    conf.face_alpha=0.5;
end
if ~isfield(conf, 'mean_colour')
    conf.mean_colour='r';
end
if ~isfield(conf, 'median_colour')
    conf.median_colour='k';
end
if ~isfield(conf, 'support_range')
    conf.support_range = conf.y_range;
end
    
% Make sure there is a face colour for each violin
if size(conf.face_color,1)==1
    conf.face_color = repmat(conf.face_color, num_violins, 1);
end

end
    
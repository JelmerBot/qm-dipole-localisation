function [h] = plotBarStackGroups(stack_data, group_labels, bar_labels)
% Plot groups of stacked bars
%
% Params: 
%      stack_data is a 3D matrix (i.e., stack_data(i, j, k) => (Group, Bar, Stack))
%      group_labels is a cell-array with one label for each group
%      bar_labels is a cell-array with one label for each bar in a group
%
% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu)
% - Updated by Jelmer Bot (2019)

%% Validate input

num_groups = size(stack_data, 1);
num_bars_per_group = size(stack_data, 2);
% numStacks = size(stackData, 3);

if length(bar_labels) ~= num_bars_per_group
    error('Then number of bar_labels does not match the number of bars per group');
end

if length(group_labels) ~= num_groups
   error('The number of group_labels does not match the number of groups'); 
end

%% Configurable stuff
% Fraction of 1. If 1, then we have all groups touching
max_group_width = 0.8;
% Fraction of 1. If 1, then we have all bars touching
max_bar_width = 0.8;
% Width of a character is cm. (Measure or estimate based on font settings)
% char_width = 0.15; %cm

%% Internal values
group_offset = 1:num_groups;
bar_width = max_group_width / num_bars_per_group;

% Keep track of the positions
group_draw_positions = zeros(num_bars_per_group, num_groups);

%% Plot the bars
hold off
for i=1:num_bars_per_group
    Y = squeeze(stack_data(:,i,:));
    
    % Center the bars:
    internal_pos_count = i - ((num_bars_per_group+1) / 2);
    
    % Offset the group draw positions:
    group_draw_positions(i, :) = (internal_pos_count) * bar_width + group_offset;
    set(gca,'ColorOrderIndex',1)
    h(i,:) = bar(Y, 'stacked');
    hold on
    set(h(i,:),'BarWidth', max_bar_width * bar_width);
    set(h(i,:),'XData',group_draw_positions(i,:));
end
hold off;

%% Create the labels and ticks
half_group_size = max_group_width / 2;
spacing = 1-max_group_width;
xlim([1 - half_group_size - spacing, num_groups + half_group_size + spacing]); 
ax = gca;
ax.TickDir = 'out';                % Stylistic choice
ax.TickLabelInterpreter = 'Latex'; % Required for multi-line tick-labels

make_label = @(l1, l2) ['\begin{tabular}{c} ' l1 ' \\ ' l2 ' \end{tabular}'];

if mod(num_bars_per_group, 2) == 1 % isodd 
    % The center tick, which gets the group-label is contained in
    % groupDrawPositions. So we only need to prepare the labels
    half = floor(num_bars_per_group / 2);
    center = half + 1;
    labels = cell(1, num_groups * num_bars_per_group);
    
    for idx = 1:length(labels)
        bar_idx = mod(idx - 1, num_bars_per_group) + 1;
        group_idx = floor((idx - 1) / num_bars_per_group) + 1;
        
        if bar_idx == center
            labels(idx) = {make_label(bar_labels{bar_idx}, group_labels{group_idx})};
        else
            labels(idx) = {make_label(bar_labels{bar_idx}, ' ')};
        end
    end
    
    ax.XTick = group_draw_positions(:);
    ax.XTickLabels = labels;
else
    % The group-tick is placed between bars, so add it manually
    num_ticks = num_groups * (num_bars_per_group + 1);
    positions = zeros(1, num_ticks);
    half = floor(num_bars_per_group / 2);
    center = half + 1;
    group_tick_idx = center:(num_bars_per_group + 1):num_ticks;
    positions(group_tick_idx) = 1:num_groups;
    positions(positions == 0) = group_draw_positions(:);
    
    % Prepare the labels
    labels = cell(1, num_groups * (num_bars_per_group + 1));
    for idx = 1:length(labels)
        bar_idx = mod(idx - 1, (num_bars_per_group + 1)) + 1;
        group_idx = floor((idx - 1) / (num_bars_per_group + 1)) + 1;
        
        if bar_idx == center
            % The center tick for the group label
            labels(idx) = {make_label(' ', group_labels{group_idx})};
        elseif bar_idx < center
            % Base case
            labels(idx) = {make_label(bar_labels{bar_idx}, ' ')};
        else
            % We insterted a tick manually, so correct the bar_index
            labels(idx) = {make_label(bar_labels{bar_idx - 1}, ' ')};
        end
    end
    
    ax.XTick = positions;
    ax.XTickLabels = labels;
end

end 

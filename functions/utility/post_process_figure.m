function post_process_figure(factor, ratio, axis_size, margin)
%post_process_figure
% A function that resizes figures based on a hardcoded maximum width
%   factor -    The factor of the width to use
%   ratio -     The aspect ratio to create (height / width)
%   axis_size - The size of the axis [y, x] in cm
%               The function will check their position and correct
%               accordingly. Just measure them on screen!
%   [margin] -  The margin [width, height] to be added opposite of the axis
%               in cm. If there are no axis or the axis are contained in
%               the image, then the margin is applied to both sides.

narginchk(3, 4)

if nargin < 4
    margin = [0 0]; % cm
end

fig = gcf;
try
    set(fig.Children, 'TickLabelInterpreter', 'latex')
catch e
    warning('Could not set TickLabelInterpreter on children')
end
set(fig, 'InvertHardcopy', 'off');
set(fig, 'PaperPositionMode', 'auto')
set(fig, 'Color', 'w')
set(fig.Children,'Box','off')

% Make sure the graphics are update before continuing
drawnow 

% Output text_width
text_width = 15.59778; % cm
% Compute the figure size
fig_width = factor * text_width;
fig_height = ratio * fig_width;

% Warn about outside legends
legend_handle = findobj(gcf, 'Type', 'Legend');
if ~isempty(legend_handle) && (...
     strcmp(legend_handle.Location, 'northoutside') ||...
     strcmp(legend_handle.Location, 'eastoutside') ||...
     strcmp(legend_handle.Location, 'westoutside') ||...
     strcmp(legend_handle.Location, 'bestoutside'))
    warning('Legend may be cutt of by the resize')
end

% Update the figure size
set(gca, 'Units', 'normalized');
fig_pos(1:2) = [5 5];
fig_pos(3) = fig_width;
fig_pos(4) = fig_height;
set(gcf, 'Units', 'centimeters', 'Position', fig_pos)
% Set paper size to be figure size for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', fig_pos);

% Update the axis
ax_visibility = get(gca, 'Visible');

% Check for subplots
h = get(gcf, 'Children');
indices = arrayfun(@(str) strcmp(str, ''), get(h, 'tag'));
if isempty(indices) || sum(indices) == 1 % There are no subplots
    if strcmp(ax_visibility, 'off')
        % Axis are off so we can scale the content to the entire screen
        set(gca, 'Units', 'centimeters');
        ax_pos = [margin(1), margin(2), fig_pos(3)-2*margin(1), fig_pos(4)-1*margin(2)];
        set(gca, 'Position', ax_pos, 'Units', 'normalized')
    else
        % Axis are on so we have to make sure they stay on screen.
        set(gca, 'Units', 'centimeters');
        ax_pos = get(gca, 'Position');   

        x_loc = get(gca, 'XAxisLocation');
        if strcmp(x_loc, 'bottom')
            ax_pos(2) = axis_size(2);
            ax_pos(4) = fig_pos(4) - axis_size(2) - margin(2);
        elseif strcmp(x_loc, 'top')
            ax_pos(2) = margin(2);
            ax_pos(4) = fig_pos(4) - axis_size(2) - margin(2);
        else
            % The axis is at origin
            % So the position depends on the scale of the axis
            % We need to check how much of the axis is outside of the content
            y_lim = get(gca, 'YLim');
            x_scale = get(gca, 'XScale');
            y_dir = get(gca, 'YDir');
            if ~strcmp(x_scale, 'linear')
               warning('non-linear scales not supported') 
            end

            % Assume the axis is within the content
            height = fig_pos(4);
            ax_pos(2) = margin(2);
            ax_pos(4) = height - margin(2);

            percentage_along_axis = position_in(0, y_lim(1), y_lim(2));

            height_bottom = height * percentage_along_axis;
            height_top = height * (1 - percentage_along_axis);
            if strcmp(y_dir, 'reverse')
               [height_bottom, height_top] = swap(height_bottom, height_top);
            end

            if height_bottom < axis_size(2)
               ax_pos(2) = axis_size(2) - height_bottom;
               ax_pos(4) = height - ax_pos(2) - margin(2);
            end
            if height_top < axis_size(2)
               ax_pos(2) = margin(2);
               ax_pos(4) = height - axis_size(2) + height_top - margin(2);
            end
        end

        y_loc = get(gca, 'YAxisLocation');
        if strcmp(y_loc, 'left')
            ax_pos(1) = axis_size(1);
            ax_pos(3) = fig_pos(3) - axis_size(1) - margin(1);
        elseif strcmp(y_loc, 'right')
            ax_pos(1) = margin(1);
            ax_pos(3) = fig_pos(3) - axis_size(1) - margin(1);
        else
            % The axis is at origin
            % So the position depends on the scale of the axis
            % We need to check how much of the axis is outside of the content
            x_lim = get(gca, 'XLim');
            y_scale = get(gca, 'YScale');
            x_dir = get(gca, 'XDir');
            if ~strcmp(y_scale, 'linear')
               warning('non-linear scales not supported') 
            end

            % Assume the axis is within the content
            width = fig_pos(3);
            ax_pos(1) = margin(1);
            ax_pos(3) = width - margin(1);

            percentage_along_axis = position_in(0, x_lim(1), x_lim(2));
            width_left = width * percentage_along_axis;
            width_right = width * (1 - percentage_along_axis);
            if strcmp(x_dir, 'reverse')
               [width_left, width_right] = swap(width_left, width_right); 
            end

            if width_left < axis_size(1)
               ax_pos(1) = axis_size(1) - width_left;
               ax_pos(3) = width - ax_pos(1) - margin(1);
            end
            if width_right < axis_size(1)
               ax_pos(1) = margin(1);
               ax_pos(3) = height - axis_size(1) + width_right - margin(1);
            end
        end

        set(gca, 'Position', ax_pos, 'Units', 'normalized');
    end
end

% Make sure the graphics are update before continuing
drawnow

end

function res = place_in(x, low, high)
    % checks whether x is in the range [low high]
    % -1 is smaller than low
    %  0 is within [low high]
    %  1 is larger than high
    if x < low
        res = -1;
        return;
    end
    if x > high
        res = 1;
        return
    end
    res = 0;
end

function res = position_in(x, low, high)
    % returns the position of x in the range [low high]
    res = place_in(x, low, high);
    if res ~= 0
        return
    end
    
    low_dist = x - low;
    high_dist = high - x;
    res = low_dist / (low_dist + high_dist);
end

function [y, x] = swap(x, y)
end
function plot_predicted_locations_quad(predictions, sources, in_focus_area, params)
%plot_predicted_locations_quad
%
% Syntax: plot_predicted_locations_quad(predictions, sources, in_focus_area params)
%
% Creates a connected locations plot of the actual and predicted locations.
narginchk(4, 4)
nargoutchk(0, 0)
if any(size(predictions) ~= size(sources))
    error('There should be the same number of predictions as there are sources');
end

% n_colors = 1024;
% cmax = pi;
% cmin = 0;
% index = fix((mod(sources(:,3),pi) - cmin)/(cmax-cmin)*(n_colors-1))+1;
% orig_map = redblue2(n_colors);
% orig_colors = orig_map(index,:);
indices = 1:size(predictions, 1);

x_range = params.domain.x_range;
y_range = params.domain.y_range;

% figure(2)
pause(0.5)
drawnow
hold off
drawnow
pause(0.5)

for sample = indices
    if numel(sample) > 1
        error('wrong loop variable')
    end
    plot(...
        [sources(sample, 1); predictions(sample, 1)],...
        [sources(sample, 2); predictions(sample, 2)],...
        '-', 'color', [0 0 0, 0.2]);%[orig_colors(sample,:) 0.1]);
    hold on
end 
hold on
p = plot(predictions(:, 1), predictions(:, 2),...
     '.', 'markersize', 4);

drawnow
markers = p.MarkerHandle;
markers.FaceColorData = uint8(255 * [0 0 0 0.7])';
markers.EdgeColorData = uint8(255 * [0 0 0 0.7])';
set(gca, 'ColorOrderIndex', 1)
plot(sources(in_focus_area, 1), sources(in_focus_area, 2), 'o',...
     'markersize', 3)
set(gca, 'ColorOrderIndex', 2)
plot(sources(~in_focus_area, 1), sources(~in_focus_area, 2), 'o',...
     'markersize', 3)
 
mk = plot(params.sensors.locations.x, 0, 'ko', 'markerfacecolor', [0 0 0]);
plot(x_range, repmat(y_range(1), 2, 1), '-', 'color', [0 0 0 0.5])
plot(x_range, repmat(y_range(2), 2, 1), '-', 'color', [0 0 0 0.5])
plot(repmat(x_range(1), 2, 1), y_range, '-', 'color', [0 0 0 0.5])
plot(repmat(x_range(2), 2, 1), y_range, '-', 'color', [0 0 0 0.5])

axis image
xlim(x_range)
ylim([-0.025 y_range(end)])

% Ticks
set(gca, 'TickDir', 'out')
% xticks(-2:0.1:2);
% xticklabels({'' '-0.4' '' '-0.2' '' '0' '' '0.2' '' '0.4' ''}); 
yticks(0:0.1:1)
yticklabels({'0' '' '0.2' '' '0.4' '' '0.6' '' '0.8' '' '1'}); 

% Labels
xlabel('$x$~(m)')
ylabel('$y$~(m)')

% colormap(redblue2(1024))
% c = colorbar;
% c.Label.Interpreter = 'Latex';
% c.Label.String = '$\varphi$~(rad)';
drawnow
post_process_figure(1, 0.45, [1.15 1], [0.3 0.1])
curunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
cursize = get(gca, 'Position');
set(gca, 'Units', curunits);
pt_sz = diff(xlim) / cursize(3);
sensor_sz = 8e-3 / pt_sz;
set(mk, 'MarkerSize', sensor_sz)
% addpi = @(c) ['$', c{1}, '\pi$'];
% c.TickLabels(2:end) = arrayfun(addpi, c.TickLabels(2:end), 'UniformOutput', false);
% 
% for sample = indices
%     quiver(sources(sample, 1), sources(sample, 2),...
%            cos(sources(sample, 3)), sin(sources(sample, 3)),...
%            10 * pt_sz, 'k-', 'color', orig_colors(sample, :));
% end

end
function plot_orientation_block(a, y, data, params, conf)

azi_in_degrees = -a / pi * 180 + 90;

% process data
color_index = zeros(size(data, 1), 1);
lvls = color_index;
for idx = 1:size(data, 1)
    lvl = find(data(idx,:) <= conf.thresholds, 1, 'first');
    if isempty(lvl)
        lvl = 6;
    end
    color_index(idx) = lvl - 1; % tresholds starts at 0!
    lvls(idx) = lvl;
end
color_index = reshape(color_index, [length(unique(y)), length(unique(a))]);

%compute colormap
my_map = flip(conf.colormap);
cmax = conf.thresholds(end-1);
cmin = 0;
index = fix((conf.thresholds(1:end-1) - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
map = squeeze(ind2rgb(index', my_map));

% make the plot
figure
polarPcolor(y', azi_in_degrees', color_index,...
    'colormap', map,...
    'NSpokes', conf.spokes,...
    'NCircles', conf.circles,...
    'colBar', 0);
post_process_figure(0.4, 0.65, [0 0], [0.4 -0.3])
text(-0.15, 1.26, '$\varphi$ (rad)')
text(0.3, -0.1', '$y$ (m)')
caxis([1 size(map, 1)])

drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
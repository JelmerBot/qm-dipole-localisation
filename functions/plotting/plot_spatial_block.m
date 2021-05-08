function plot_spatial_block(x, y, data, params, conf)

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
color_index = reshape(color_index, [length(unique(y)), length(unique(x))]);

% compute colormap
my_map = flip(conf.colormap);
cmax = conf.thresholds(end-1);
cmin = 0;
index = fix((conf.thresholds(1:end-1) - cmin)/(cmax-cmin)*(length(my_map)-1))+1;
map = squeeze(ind2rgb(index', my_map));

% make the plot
figure
hold off
imagesc(x, y, color_index)
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size / 2 max(y) + conf.cell_size / 2])
colormap(map)
caxis([1 size(map,1)])
daspect([1 1 1])
xlabel('$x$ (m)')
ylabel('$y$ (m)')

post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
set(gca, 'tickdir', 'out')
drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
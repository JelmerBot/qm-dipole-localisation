function plot_spatial_heatmap(x, y, data, params, conf)

data = reshape(data, [length(y), length(x)]);

figure
hold off
imagesc(x, y, data)
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size / 2 max(y) + conf.cell_size / 2])
colormap(conf.colormap)
caxis(conf.caxis)
daspect([1 1 1])
xlabel('$x$ (m)')
ylabel('$y$ (m)')

% c = colorbar;
% c.Label.Interpreter = 'Latex';
% c.Label.String = 'error score';
% c.Label.FontSize = 8;

post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
set(gca, 'tickdir', 'out')
drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
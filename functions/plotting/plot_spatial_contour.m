function plot_spatial_contour(x, y, data, params, conf)

data = reshape(data, [length(y), length(x)]);
[X, Y] = meshgrid(x, y);

figure
hold off
contourf(X, Y, data, conf.thresholds)
set(gca, 'ydir', 'normal')
hold on
plot(params.sensors.locations.x, params.sensors.locations.y, 'k.');
ylim([0-conf.cell_size/2 max(y)])

xlabel('$x$ (m)')
ylabel('$y$ (m)')

colormap(conf.colormap)
caxis([0 conf.thresholds(end-1)])
% c = colorbar;
% c.Label.Interpreter = 'Latex';
% c.Label.String = 'median Error';
% c.Label.FontSize = 8;
% c.Ticks = conf.thresholds;
% c.TickLength = 0.13;

daspect([1 1 1])

post_process_figure(0.4, 0.63, [0.95 0.75], [0.3 0])
set(gca, 'tickdir', 'out')
drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end


function plot_orientation_contour(a, y, data, params, conf)

azi_in_degrees = -a / pi * 180 + 90;
data = reshape(data, [length(y), length(a)]);

figure
polarPcontour(y', azi_in_degrees', data, conf.thresholds,...
    'colormap', conf.colormap,...
    'NSpokes', conf.spokes,...
    'NCircles', conf.circles,...
    'colBar', 0);
caxis(conf.caxis)
post_process_figure(0.4, 0.65, [0 0], [0.4 -0.3])
text(-0.15, 1.26, '$\varphi$ (rad)')
text(0.3, -0.1', '$y$ (m)')

drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
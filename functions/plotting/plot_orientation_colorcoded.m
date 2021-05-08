function plot_orientation_colorcoded(a, y, data, params, conf)

azi_in_degrees = -a / pi * 180 + 90;
data = reshape(data, [length(y), length(a)]);

figure
polarPcolor(y', azi_in_degrees', data,...
    'colormap', conf.colormap,...
    'NSpokes', conf.spokes,...
    'NCircles', conf.circles,...
    'colBar', 0);
post_process_figure(0.24, 0.7, [0 0], [0.4 -0.3])
text(-0.15, 1.46, '$\varphi$ (rad)')
text(0.3, -0.15, '$y$ (m)')
caxis([1 size(conf.colormap, 1)])

drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
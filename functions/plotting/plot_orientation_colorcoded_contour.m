function plot_orientation_colorcoded_contour(a, y, data, params, conf)

azi_in_degrees = -a / pi * 180 + 90;
data = reshape(data, [length(y), length(a)]);

data(data == 6) = 30;
data(data == 12) = 30;
data(data == 18) = 30;
data(data == 24) = 30;
data(data > 30) = 30;

indices = 1:size(conf.colormap, 1);
indices([6 12 18 24 31:end]) = [];

figure
polarPcontour(y', azi_in_degrees', data, indices,...
    'colormap', conf.colormap,...
    'NSpokes', conf.spokes,...
    'NCircles', conf.circles,...
    'colBar', 0);
post_process_figure(0.4, 0.65, [0 0], [0.4 -0.3])
text(-0.15, 1.26, '$\varphi$ (rad)')
text(0.3, -0.1', '$y$ (m)')
caxis([1 size(conf.colormap, 1)])

drawnow
pause(0.2)
print(fullfile(conf.name), '-djpeg', '-r600')
disp(conf.name)
close

end
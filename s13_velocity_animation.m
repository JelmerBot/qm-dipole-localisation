%% Introduction
% Generates a 3D velocity animation. This animation was used during the
% design of the QM method.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Configure animation

seconds = 7;
fps = 30;
step = 1/fps;
steps = seconds * fps;
periods = 2;
animation_filename = 'images/velocity_animation/rotating_velocity.gif';
animation_begin = 3;
animation_end = 4;
animation_start_pos = [0 90];
animation_end_pos = [-30 30];

%% Compute animation viewpoint

% Compute viewing angle
animation_start_step = animation_begin * fps;
animation_end_step = animation_end * fps;
n_animation_steps = animation_end_step - animation_start_step;
panning = interp1([animation_start_step; animation_end_step],...
                  [animation_start_pos; animation_end_pos],...
                  animation_start_step:animation_end_step);
views = [repmat(animation_start_pos, animation_start_step - 1, 1);...
         panning;...
         repmat(animation_end_pos, steps - animation_end_step, 1)];

%% Compute wavelets and envelope

N=1024;                                             % #of points used for wavelet 
b=0;                                                % Source location along x-axis        
rhorange=5;                                         % Range of rho values (from - to +)
C=1;                                                % Overall normalization constant
lin=linspace(1,2*N+1,2*N+1);                        % Creates an array (linear increase without offset)
x=lin-N-1 ;                                         % Symmetric x-vector
d=N/rhorange;                                       % Source distance to array
rho=(x-b)/d;                                        % Normalization of x with respect to distance
x_plot = linspace(-rhorange, rhorange, length(x));

% The orientations of the source
phis = linspace(0, periods*2*pi, steps);
phis = mod(phis, 2*pi);

% wavelet values
denom=(1+rho.^2).^(5/2);                        %Denominator of the wavelets
psi_e=(-1+2*rho.^2)./denom;                     %Even wavelet
psi_o=(-3*rho)./denom;                          %Odd wavelet
psi_n=(2-rho.^2)./denom;                        %Navelet
% envenope values
env_vx=sqrt(psi_e.^2+psi_o.^2);                 %Envelope of vx
env_vy=sqrt(psi_o.^2+psi_n.^2);                 %Envelope of vy
% easy 3d plotting value
zero = zeros(size(x));

%% Render animation

figure(5)
colors = get(gca, 'ColorOrder');
post_process_figure(1, 0.4, [0.5 0.5], [0 0])
first = true;
for t = 1:steps
    tic;
    phi=phis(t);
    
    subplot(131)
    hold off
    %plot envelope
%     plot3(x_plot, env_vx, zero, 'g-')
%     hold on
%     plot3(x_plot, -env_vx, zero, 'g-')
%     plot3(x_plot, zero, env_vx, 'g-')
%     plot3(x_plot, zero, -env_vx, 'g-')
    % plot velocity
    p = [x_plot', psi_e', psi_o'];
    p_rot = arrayfun(@(idx) rotx(phi*180/pi) * p(idx,:)', 1:size(p,1), 'UniformOutput', false);
    p_rot = [p_rot{:}]';
    set(gca, 'ColorOrderIndex', 1)
    plot3(p_rot(:,1), p_rot(:,2), p_rot(:,3), '-');
    grid on
    zlim([-2 2])
    ylim([-2 2])
    view(views(t,1), views(t,2))
    xlabel('$x$')
    ylabel('$\varphi_e$')
    zlabel('$\varphi_o$')
    daspect([1 1 1])
    title('$v_x$ ')
    
    subplot(132)
    hold off
%     %plot envelope
%     plot3(x_plot, env_vy, zero, 'g-')
%     hold on
%     plot3(x_plot, -env_vy, zero, 'g-')
%     plot3(x_plot, zero, env_vy, 'g-')
%     plot3(x_plot, zero, -env_vy, 'g-')
    % plot velocity
    p = [x_plot', psi_o', psi_n'];
    p_rot = arrayfun(@(idx) rotx(phi*180/pi) * p(idx,:)', 1:size(p,1), 'UniformOutput', false);
    p_rot = [p_rot{:}]';
    plot3(p_rot(:,1), p_rot(:,2), p_rot(:,3), '-', 'color', colors(2,:));
    grid on
    zlim([-2 2])
    ylim([-2 2])
    view(views(t,1), views(t,2))
    xlabel('$x$')
    ylabel('$\varphi_o$')
    zlabel('$\varphi_n$')
    daspect([1 1 1])
    title('$v_y$')
    
    subplot(133)
    hold off
    polarplot([phi, phi], [0 1], 'color', colors(1,:));
    hold on
    rticks([])
    title('$\varphi$ (rad)')
    ax = gca;
    ax.ThetaAxis.TickLabelInterpreter = 'latex';
    thetaticklabels({'$0$', '', '', '$\frac{1}{2}\pi$', '' '', '$\pi$', '', '', '$1\frac{1}{2}\pi$', '' ,''})
    fig = gcf;
    polarfig = fig.Children(1);
    polarfig.Position = polarfig.Position .* [1 1 0.8 1];
    polarfig.Position(1) = polarfig.Position(1) + 0.07;
    vyfig = fig.Children(2);
    vyfig.Position(1) = vyfig.Position(1) + 0.05;
    
    drawnow
    if first
        gif(animation_filename, 'Loop', 1, 'DelayTime', step, 'frame', gcf)
        first = false;
    else
        gif()
    end
%     pause(step - toc)
end

%% y-z over lateral-line

animation_filename = 'images/velocity_animation/wavelet_contributions_over_all.gif';
seconds = 7;
fps = 30;
step = 1/fps;
steps = seconds * fps;
phi = 0.6*pi;
begin_pos = 0.3;
end_pos = 0.7;
n_indices = length(psi_e);
l_y = 2;
l_x = l_y;

indices = round(linspace(begin_pos * n_indices, end_pos * n_indices, steps));

figure(6)
post_process_figure(1, 0.4, [0 0], [0 0])
first = true;
for idx = indices
    tic
    colors = get(gca, 'ColorOrder');
    subplot(131)
    % Compute contribution of both wavelets
    [psi_e_x, psi_e_y] = pol2cart(phi, cos(phi)*psi_e(idx));
    [psi_o_x, psi_o_y] = pol2cart(phi, sin(phi)*psi_o(idx));
    x = psi_o_x + psi_e_x;
    y = psi_o_y + psi_e_y;
    
    % plot velocity vector
    hold off
    plot([0 x], [0, y], '-', 'color', colors(1,:));
    hold on

    % Plot wavelet axis
    plot([0 psi_e(idx)], [0 0], 'k-')
    plot([0 0], [0 psi_o(idx)], 'k-')

    % Plot envelope
    plot([0 psi_e(idx)], [0, psi_o(idx)], '-', 'color', colors(3,:))
    plot([psi_e(idx) psi_e(idx)], [0 psi_o(idx)], 'k:')
    plot([0 psi_e(idx)], [psi_o(idx) psi_o(idx)], 'k:')

    % Show wavelet contributions
    plot([psi_o_x 0], [psi_o_y psi_o(idx)], 'k--')
    plot([psi_e_x psi_e(idx)], [psi_e_y 0], 'k--')

    % Show envelope contribution
    plot([psi_e(idx) x], [psi_o(idx) y], 'k--')
    daspect([1 1 1])
    axis([-l_x l_x -l_x l_x])
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')
    title('$v_x$')
    xlabel('$\psi_e$')
    ylabel('$\psi_o$')

    subplot(132)
    hold off
    [psi_o_x, psi_o_y] = pol2cart(phi, cos(phi)*psi_o(idx));
    [psi_n_x, psi_n_y] = pol2cart(phi, sin(phi)*psi_n(idx));
    x = psi_o_x + psi_n_x;
    y = psi_o_y + psi_n_y;

    % plot velocity vector
    plot([0 x], [0, y], '-', 'color', colors(2,:));
    hold on

    % Plot wavelet axis
    plot([0 psi_o(idx)], [0 0], 'k-')
    plot([0 0], [0 psi_n(idx)], 'k-')

    % Plot envelope
    plot([0 psi_o(idx)], [0, psi_n(idx)], '-', 'color', colors(3,:))
    plot([psi_o(idx) psi_o(idx)], [0 psi_n(idx)], 'k:')
    plot([0 psi_o(idx)], [psi_n(idx) psi_n(idx)], 'k:')

    % Show wavelet contributions
    plot([psi_n_x 0], [psi_n_y psi_n(idx)], 'k--')
    plot([psi_o_x psi_o(idx)], [psi_o_y 0], 'k--')

    % Show envelope contribution
    plot([psi_o(idx) x], [psi_n(idx) y], 'k--') 
    daspect([1 1 1])
    axis([-l_y l_y -l_y l_y]);
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')
    title('$v_y$')
    xlabel('$\psi_o$')
    ylabel('$\psi_n$')
    
    % Orientation plot
    subplot(233)
    hold off
    polarplot([phi, phi], [0 1]);
    hold on
    rticks([])
    
    ax = gca;
    ax.ThetaAxis.TickLabelInterpreter = 'latex';
    ax.ThetaAxis.Label.String = '$\varphi$ (rad)';
    thetaticklabels({'$0$', '', '', '$\frac{1}{2}\pi$', '' '', '$\pi$', '', '', '$1\frac{1}{2}\pi$', '' ,''})
    
    % Indication of position
    subplot(236)
    hold off
%     plot(x_plot(indices), 0, 'k.')
    plot(x_plot, psi_e, 'k-')
    hold on
    plot(x_plot, psi_o, 'k-')
    plot(x_plot, psi_n, 'k-')
    plot(x_plot(idx), 0, '.', 'MarkerSize', 15)
    hold on
    daspect([1 1 1])
    xlabel('$x$')
%     xlim([min(x_plot), max(x_plot)])
%     ylim([-0.01 0.01])
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')

    
    drawnow
    if first
        gif(animation_filename, 'Loop', 1, 'DelayTime', step, 'frame', gcf)
        first = false;
    else
        gif()
    end
%     pause(step-toc)
end


%% y-z over orientation

animation_filename = 'images/velocity_animation/wavelet_contributions_over_phi.gif';
seconds = 7;
fps = 30;
step = 1/fps;
steps = seconds * fps;
begin_pos = 0.55;
idx = round(length(psi_e) * begin_pos);
l_y = 1.2;
l_x = l_y;
periods = 1;
phis = linspace(0, periods*2*pi, steps);
phis = mod(phis, 2*pi);

figure(6)
post_process_figure(1, 0.4, [0 0], [0 0])
first = true;
for phi = phis
    tic
    colors = get(gca, 'ColorOrder');
    subplot(131)
    % Compute contribution of both wavelets
    [psi_e_x, psi_e_y] = pol2cart(phi, cos(phi)*psi_e(idx));
    [psi_o_x, psi_o_y] = pol2cart(phi, sin(phi)*psi_o(idx));
    x = psi_e_x + psi_o_x;
    y = psi_e_y + psi_o_y;
    
    % plot velocity vector
    hold off
    plot([0 x], [0, y], '-', 'color', colors(1,:));
    hold on

    % Plot wavelet axis
    plot([0 psi_e(idx)], [0 0], 'k-')
    plot([0 0], [0 psi_o(idx)], 'k-')

    % Plot envelope
    plot([0 psi_e(idx)], [0, psi_o(idx)], '-', 'color', colors(3,:))
    plot([psi_e(idx) psi_e(idx)], [0 psi_o(idx)], 'k:')
    plot([0 psi_e(idx)], [psi_o(idx) psi_o(idx)], 'k:')

    % Show wavelet contributions
    plot([psi_o_x 0], [psi_o_y psi_o(idx)], 'k--')
    plot([psi_e_x psi_e(idx)], [psi_e_y 0], 'k--')

    % Show envelope contribution
    plot([psi_e(idx) x], [psi_o(idx) y], 'k--')
    daspect([1 1 1])
    axis([-l_x l_x -l_x l_x])
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')
    title('$v_x$')
    xlabel('$\psi_e$')
    ylabel('$\psi_o$')

    subplot(132)
    hold off
    [psi_o_x, psi_o_y] = pol2cart(phi, cos(phi)*psi_o(idx));
    [psi_n_x, psi_n_y] = pol2cart(phi, sin(phi)*psi_n(idx));
    x = psi_o_x + psi_n_x;
    y = psi_o_y + psi_n_y;

    % plot velocity vector
    plot([0 x], [0, y], '-', 'color', colors(2,:));
    hold on

    % Plot wavelet axis
    plot([0 psi_o(idx)], [0 0], 'k-')
    plot([0 0], [0 psi_n(idx)], 'k-')

    % Plot envelope
    plot([0 psi_o(idx)], [0, psi_n(idx)], '-', 'color', colors(3,:))
    plot([psi_o(idx) psi_o(idx)], [0 psi_n(idx)], 'k:')
    plot([0 psi_o(idx)], [psi_n(idx) psi_n(idx)], 'k:')

    % Show wavelet contributions
    plot([psi_n_x 0], [psi_n_y psi_n(idx)], 'k--')
    plot([psi_o_x psi_o(idx)], [psi_o_y 0], 'k--')

    % Show envelope contribution
    plot([psi_o(idx) x], [psi_n(idx) y], 'k--') 
    daspect([1 1 1])
    axis([-l_y l_y -l_y l_y]);
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')
    title('$v_y$')
    xlabel('$\psi_o$')
    ylabel('$\psi_n$')
    
    % Orientation plot
    subplot(233)
    hold off
    polarplot([phi, phi], [0 1]);
    hold on
    rticks([])
    
    ax = gca;
    ax.ThetaAxis.TickLabelInterpreter = 'latex';
    ax.ThetaAxis.Label.String = '$\varphi$ (rad)';
    thetaticklabels({'$0$', '', '', '$\frac{1}{2}\pi$', '' '', '$\pi$', '', '', '$1\frac{1}{2}\pi$', '' ,''})
    
    % Indication of position
    subplot(236)
    hold off
    plot(x_plot, psi_e, 'k-')
    hold on
    plot(x_plot, psi_o, 'k-')
    plot(x_plot, psi_n, 'k-')
    plot(x_plot(idx), 0, '.', 'MarkerSize', 15)
    hold on
    daspect([1 1 1])
    xlabel('$x$')
%     xlim([min(x_plot), max(x_plot)])
%     ylim([-0.01 0.01])
    set(gca, 'tickdir', 'out')
    set(gca,'Box','off')

    
    drawnow
    if first
        gif(animation_filename, 'Loop', 1, 'DelayTime', step, 'frame', gcf)
        first = false;
    else
        gif()
    end
%     pause(step-toc)
end

%% 3D Viewpoint images

phi = 0;
output_folder = 'images/velocity_animation';

% vx
p = [x_plot', psi_e', psi_o'];
p_rot = arrayfun(@(idx) rotx(phi*180/pi) * p(idx,:)', 1:size(p,1), 'UniformOutput', false);
p_rot = [p_rot{:}]';

figure
plot3(p_rot(:,1), p_rot(:,2), p_rot(:,3), '-');
grid on
zlim([-2 2])
ylim([-2 2])
xlabel('$\rho$')
ylabel('$\psi_e$')
zlabel('$\psi_o$')
daspect([1 1 1])
set(gca, 'tickdir', 'out')

view([1 0 0])
post_process_figure(0.3, 1, [1 1], [0.5 0.5])
print([output_folder, '/velocity_x_yz.jpg'], '-djpeg', '-r600')

xlim([-2 2])
view(2)
post_process_figure(0.3, 1, [1 1], [0.5 0.5])
print([output_folder, '/velocity_x_xy.jpg'], '-djpeg', '-r600')

view(3)
post_process_figure(0.3, 1, [0.8 0.8], [0.5 0.1])
print([output_folder, '/velocity_x_3d.jpg'], '-djpeg', '-r600')
close

% vy
p = [x_plot', psi_o', psi_n'];
p_rot = arrayfun(@(idx) rotx(phi*180/pi) * p(idx,:)', 1:size(p,1), 'UniformOutput', false);
p_rot = [p_rot{:}]';

figure

colors = get(gca, 'ColorOrder');
plot3(p_rot(:,1), p_rot(:,2), p_rot(:,3), '-', 'color', colors(2,:));
grid on
zlim([-2 2])
ylim([-2 2])
xlabel('$\rho$')
ylabel('$\varphi_o$')
zlabel('$\varphi_n$')
daspect([1 1 1])
set(gca, 'tickdir', 'out')

view([1 0 0])
post_process_figure(0.3, 1, [1 1], [0.5 0.5])
print([output_folder, '/velocity_y_yz.jpg'], '-djpeg', '-r600')

xlim([-2 2])
view(2)
post_process_figure(0.3, 1, [1 1], [0.5 0.5])
print([output_folder, '/velocity_y_xy.jpg'], '-djpeg', '-r600')

view(3)
post_process_figure(0.3, 1, [0.8 0.8], [0.5 0.1])
print([output_folder, '/velocity_y_3d.jpg'], '-djpeg', '-r600')
close

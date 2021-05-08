function plot_cwt(X, Y, Z, x_even, y_nav, x_odd, y_odd, p, p_hat, params)

c_x = params.cwt.c_x;
c_y = params.cwt.c_y;

subplot(5,2,1)
plot_cwt_scatter(X, Y, x_even, p, p_hat, params)
colormap(gca, redblue(2014));
title(gca, '$x_{even}$')
subplot(5,2,2)
plot_cwt_scatter(X, Y, y_nav, p, p_hat, params)
colormap(gca, redblue(2014));
title(gca, '$y_{nav}$')

subplot(5,2,3)
plot_cwt_scatter(X, Y, x_odd, p, p_hat, params)
colormap(gca, redblue(2014));
title(gca, '$x_{odd}$')
subplot(5,2,4)
plot_cwt_scatter(X, Y, y_odd, p, p_hat, params)
colormap(gca, redblue(2014));
title(gca, '$y_{odd}$')

subplot(5,2,5)
plot_cwt_scatter(X, Y, sqrt(x_even.^2 + x_odd.^2), p, p_hat, params)
title(gca, '$x_{mag}$')
colormap(gca, magma(2014));
subplot(5,2,6)
plot_cwt_scatter(X, Y, sqrt(y_nav.^2 + y_odd.^2), p, p_hat, params)
colormap(gca, magma(2014));
title(gca, '$y_{mag}$')

subplot(5,2,7)
plot_cwt_scatter(X, Y, sqrt(x_even.^2 + c_x * x_odd.^2), p, p_hat, params)
colormap(gca, magma(2014));
title(gca, '$x_{mag}$ corrected')
subplot(5,2,8)
plot_cwt_scatter(X, Y, sqrt(c_y * y_nav.^2 + y_odd.^2), p, p_hat, params)
colormap(gca, magma(2014));
title(gca, '$x_{mag}$ corrected')


subplot(5,2,9)
plot_cwt_scatter(X, Y, sqrt(x_even.^2 + c_x * x_odd.^2 + c_y * y_nav.^2 + y_odd.^2), p, p_hat, params)
colormap(gca, magma(2014));
title(gca, '$x+y$ corrected')
subplot(5,2,10)
plot_cwt_scatter(X, Y, sqrt(x_even.^2 + x_odd.^2 + y_nav.^2 + y_odd.^2), p, p_hat, params)
colormap(gca, magma(2014));
title(gca, '$x+y$')

end

function plot_cwt_scatter(x, y, c, p, p_hat, params)
    hold off
    scatter(x, y, 5, c, 'filled');
    hold on
    if ~isempty(p)
        plot(p(1), p(2), 'kx');
    end
%     plot(p_hat(1), p_hat(2), 'go');
%     [~, i] = max(c);
    
    p_hat = fit_guass(x, y, c, params);
    plot(p_hat(1), p_hat(2), 'go');
    
    if any(c < 0)
        c = max(abs(c(:)));
        caxis([-c c]);
    end
    xlim(params.domain.x_range);
    ylim(params.domain.y_range);
    daspect([1 1 1]);
    xlabel('$x$')
    ylabel('$y$')
end

function p_hat = fit_guass(X, Y, Z, params)

Z = Z ./ max(abs(Z));
thresholds = [params.cwt.threshold_min, params.cwt.threshold_max];
filter = Z > thresholds(1) & Z < thresholds(2);
n_sources = sum(filter(:));
if n_sources <= 0
%     warning('CWT:filterExhausted', ...
%         ['No coefficient values passed the threshold filter. '...
%          'Used the maximum value instead.'])
    [~, i] = max(X(:));
    p_hat = [X(i), Y(i)];
elseif n_sources < 7 % The number of variables that is fitted
%     warning('CWT:filterTooRestrictive', ...
%         ['Not enough points passed the threshold filter to fit a Gaussian.\n',...
%          'Used the maximum value within the filter instead.'])
    [~, i] = max(Z(filter));
    X_filter = X(filter);
    Y_filter = Y(filter);
    p_hat = [X_filter(i), Y_filter(i)];
else
    res = fmgaussfit(X(filter), Y(filter), Z(filter)', params);
    p_hat = res(5:6);
end

end
 
% function plot_cwt_heatmap(x_values, y_values, Z, p, p_hat, params)
%     [X, Y] = meshgrid(x_values, y_values);
%     imagesc(x_values, y_values, Z)
%     colormap(gca, redblue(1024))
%     set(gca, 'ydir', 'normal')
%     hold on
%     plot(p(1), p(2), 'kx')
%     plot(p_hat(1), p_hat(2), 'go');
%     [~, i] = max(Z(:));
%     plot(X(i), Y(i), 'bo');
%     c = max(abs(caxis));
%     caxis([-c c]);
%     xlim(params.domain.x_range);
%     ylim(params.domain.y_range);
%     daspect([1 1 1])
% end
% 
% function plot_cwt(x_values, y_values, W, p, p_hat, params)
% 
% figure
% subplot(3, 2, 1)
% plot_cwt_heatmap(x_values, y_values, W.x_even, p, p_hat, params);
% 
% subplot(3, 2, 2)
% plot_cwt_heatmap(x_values, y_values, W.y_nav, p, p_hat, params);
% 
% subplot(3, 2, 3)
% plot_cwt_heatmap(x_values, y_values, W.x_odd, p, p_hat, params);
% 
% subplot(3, 2, 4)
% plot_cwt_heatmap(x_values, y_values, W.y_odd, p, p_hat, params);
% 
% subplot(3, 1, 3)
% plot_cwt_heatmap(x_values, y_values, W.c_mag, p, p_hat, params);
% 
% end

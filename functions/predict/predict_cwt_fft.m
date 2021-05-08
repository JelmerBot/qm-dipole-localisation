% function [predictions, cwt] = predict_cwt_fft(velocity, wavelets, params)
% %predict_knn
% %
% % Syntax: predictions = predict_cwt_fft(velocity, wavelets, params)
% %
% % Implements cwt predictions.
% narginchk(3, 3)
% nargoutchk(1, 1)
% velocity.time = 0;
% validatevelocity(velocity)
% 
% %% Convienence variables
% sensors = params.sensors;
% n_samples = size(velocity.xt, 1);
% predictions = zeros(n_samples, 3);
% 
% % Apply input mode
% [vx, vy] = cwt_input_mode_velocity(velocity, params);
% vx = permute(vx, [1 3 2]);
% vy = permute(vy, [1 3 2]);
% wavelets = cwt_normalise_wavelets(wavelets);
% 
% %% Compute coefficient range
% M = size(wavelets.even, 2);
% N = size(vx, 2);
% step = diff(params.sensors.x_range) ./ N;
% minval = -step * (N/2 + M/2 - 1);
% maxval = step * (N/2 + M/2 - 1);
% x_values = minval:step:maxval;
% y_values = wavelets.d;
% [X, Y] = meshgrid(x_values, y_values);
% 
% %% Compute FFTS
% % Padding
% padded_size = M+N-1;
% 
% % Compute wavelet fft
% F_even = fft(wavelets.even, padded_size, 2);
% F_odd = fft(wavelets.odd, padded_size, 2);
% F_navelet = fft(wavelets.navelet, padded_size, 2);
% 
% % Compute velocity fft
% F_x = fft(vx, padded_size, 2); 
% F_y = fft(vy, padded_size, 2);
% scaling = 1 ./ sqrt(abs(sensors.y_value - wavelets.d));
% 
% %% Compute predictions
% threshold = params.cwt.threshold;
% parfor sample_idx = 1:n_samples
%     % Compute CWT coefficients
%     x_even = scaling .* ifft(F_x(sample_idx, :) .* F_even, [], 2, 'symmetric');
%     x_odd = scaling .* ifft(F_x(sample_idx, :) .* F_odd, [], 2, 'symmetric');
%     y_odd = scaling .* ifft(F_y(sample_idx, :) .* F_odd, [], 2, 'symmetric');
%     y_nav = scaling .* ifft(F_y(sample_idx, :) .* F_navelet, [], 2, 'symmetric');
%     
%     % Combine to 1 value
%     magnitude = sqrt(x_even.^2 + y_nav.^2 + x_odd.^2 + y_odd.^2);
%     
%     % Normalise and compute gaussian fit
%     Z = magnitude ./ max(magnitude(:));
%     filter = Z > threshold;
%     if sum(filter(:)) < 7 % The number of variables that is fitted
%         % warning('CWT:filterTooRestrictive', ...
%         %     ['The threshold value is too high, there are not enough points to fit a Gaussian.\n',...
%         %      'Used the location of the maximum value instead.'])
%         [~, i] = max(Z(:));
%         p_hat = [X(i), Y(i)];
%     else
%         res = fmgaussfit(X(filter), Y(filter), Z(filter), params);
%         p_hat = res(5:6);
%     end
%     
%     % Find coefficient values
%     [~, i] = min(sum(([X(:), Y(:)] - p_hat).^2, 2));
%     phi_x = atan2(x_odd(i), x_even(i)); 
%     phi_y = atan2(y_nav(i), y_odd(i));
%     phi = compute_angle_average([phi_x, phi_y]);
%     predictions(sample_idx, :) = [p_hat phi];
% end
% 
% end
% 
% function [vx, vy] = cwt_input_mode_velocity(velocity, params)
% 
% % All sensors measure both x and y
% if strcmp(params.sensors.input_mode, 'x+y')
%     vx = velocity.xt;
%     vy = velocity.yt;
%     return
% end
% 
% % All sensors measure x
% if strcmp(params.sensors.input_mode, 'x')
%     vx = velocity.xt;
%     vy = [];
%     return
% end
% 
% % All sensors measure x
% if strcmp(params.sensors.input_mode, 'y')
%     vy = velocity.yt;
%     vx = [];
%     return
% end
% 
% % Subsequent sensors measure x then y
% if strcmp(params.sensors.input_mode, 'x|y')
%     i1 = 1:2:(params.sensors.n_sensors);
%     i2 = 2:2:params.sensors.n_sensors;
%     vx = velocity.xt;
%     vy = velocity.yt;
%     vx(:, :, i2) = 0;
%     vy(:, :, i1) = 0; 
%     return
% end
% 
% end
% 
% function wavelets = cwt_normalise_wavelets(wavelets)
%     wavelets.even = wavelets.even ./ max(wavelets.even, [], 2);
%     wavelets.odd = wavelets.odd ./ max(wavelets.odd, [], 2);
%     wavelets.navelet = wavelets.navelet ./ max(wavelets.navelet, [], 2);
% end
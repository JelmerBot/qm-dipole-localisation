% function [even, x_odd, y_odd, navelet] = cwt_input_mode_wavelet(wavelets, params)
% 
% % Normalise over entire domain and select only sensor locations
% even = wavelets.even(:, wavelets.sensor_mask) ./ max(abs(wavelets.even), [], 2);
% wavelets.even = [];
% odd = wavelets.odd(:, wavelets.sensor_mask) ./ max(abs(wavelets.odd), [], 2);
% wavelets.odd = [];
% navelet = wavelets.navelet(:, wavelets.sensor_mask) ./ max(abs(wavelets.navelet), [], 2);
% wavelets.navelet = [];
% 
% % All sensors measure both x and y
% if strcmp(params.sensors.input_mode, 'x+y')
%     x_odd = odd;
%     y_odd = odd;
%     return
% end
% 
% % All sensors measure x
% if strcmp(params.sensors.input_mode, 'x')
%     x_odd = odd;
%     y_odd = [];
%     navelet = [];
%     return
% end
% 
% % All sensors measure x
% if strcmp(params.sensors.input_mode, 'y')
%     y_odd = odd;
%     x_odd = [];
%     even = [];
%     return
% end
% 
% % Subsequent sensors measure x then y
% if strcmp(params.sensors.input_mode, 'x|y')
%     i1 = 1:2:(params.sensors.n_sensors);
%     i2 = 2:2:params.sensors.n_sensors;
%     
%     even = even(:, i1);
%     x_odd = odd(:, i1);
%     y_odd = odd(:, i2);
%     navelet = navelet(:, i2);
%     return
% end
% 
% end
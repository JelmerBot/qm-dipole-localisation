function [snr_x, snr_y] = compute_snr(velocity, control, params)
%compute_snr
%
% Syntax: [snr_x, snr_y] = compute_snr(velocity, control, params)

narginchk(3, 3)
nargoutchk(2, 2)
validatevelocity(velocity)
validatevelocity(control)

if any(velocity.sources(:) ~= control.sources(:))
   error('Velocity and Control must have the same locations in the same order') 
end

velocity = compute_velocity_fft(velocity, params);
control = compute_velocity_fft(control, params);
i_source = velocity.source_freq_idx;

signal_power_xt = velocity.power_xt(:,i_source,:);
noise_power_xt = control.power_xt(:,i_source,:);
ratios_x = signal_power_xt ./ noise_power_xt;
snr_x = mag2db(ratios_x);

signal_power_yt = velocity.power_yt(:,i_source,:);
noise_power_yt = control.power_yt(:,i_source,:);
ratios_y =  signal_power_yt ./ noise_power_yt;
snr_y = mag2db(ratios_y);

end

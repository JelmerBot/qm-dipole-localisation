function wavelets = select_wavelet_indices(wavelets, indices)
%select_wavelet_indices
%
% Syntax: velocity = select_wavelet_indices(wavelets, indices)
narginchk(2, 2)
nargoutchk(1, 1)

wavelets.even = wavelets.even(indices, :);
wavelets.odd = wavelets.odd(indices, :);
wavelets.navelet = wavelets.navelet(indices, :);
wavelets.sources = wavelets.sources(indices,:);

end
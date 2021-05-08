function sources_bin = compute_bins(sources, cell_size)
%compute_bins
%
% Syntax: sources_bin = compute_bins(sources, cell_size)

narginchk(2,2)
nargoutchk(1,1)

scaling = 1e6;

scaled_x = scaling .* sources(:,1);
scaled_y = scaling .* sources(:,2);
scaled_azi = scaling .* sources(:,3);
scaled_cell_size = scaling .* cell_size;

x = apply_binning(scaled_x, scaled_cell_size, cell_size);
y = apply_binning(scaled_y, scaled_cell_size, cell_size);
azi = apply_binning(scaled_azi, pi.*scaled_cell_size, pi.*cell_size);
sources_bin = [x, y, azi];

end

function res = apply_binning(scaled_values, scaled_cell_size, cell_size)
    t1 = floor(scaled_values / scaled_cell_size);
    t2 = round(1 / scaled_cell_size * mod(scaled_values, scaled_cell_size));
    res = (t1 + t2) * cell_size;
end
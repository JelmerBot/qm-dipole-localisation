function weights = he(input_size, output_size)
%he
%
% Syntax: weights = he(input_size, output_size)
%
% He initialization
%
% Normal distribution with mean 0
% Variance 2 / input_size

sz = [input_size output_size];
multiplier = sqrt(2 / sum(input_size));
weights = multiplier .* randn(sz);

end
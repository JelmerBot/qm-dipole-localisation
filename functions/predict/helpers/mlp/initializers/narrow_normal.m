function weights = narrow_normal(input_size, output_size)
%narrow_normal
%
% Syntax: weights = narrow_normal(input_size, output_size)
%
% Narow normal initialization
%
% Gaussian distribution with mean 0
% standard deviation 0.01

sz = [input_size output_size];
weights = randn(sz) * 0.01;

end
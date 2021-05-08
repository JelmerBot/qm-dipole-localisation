function weights = glorot(input_size, output_size)
%glorot
%
% Syntax: weights = glorot(input_size, output_size)
%
% Grolot initialization
%
% Uniform distribution with mean 0
% Range between +/- sqrt(6 / (input_size + output_size))

sz = [input_size output_size];
multiplier = sqrt(6 / sum(sz));
weights = multiplier .* (2 .* rand(sz) - 1);

end
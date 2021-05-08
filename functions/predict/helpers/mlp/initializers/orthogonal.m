function weights = orthogonal(input_size, output_size)
%orthogonal
%
% Syntax: weights = he(input_size, output_size)
%
% Orthogonal initialization
%
% Creates a random orthogonal matrix

validateattributes(input_size, {'numeric'}, {'scalar'})
validateattributes(output_size, {'numeric'}, {'scalar'})

if output_size >= input_size
   weights = makeOrthogonalWeights(output_size, input_size)';
else
   weights = makeOrthogonalWeights(input_size, output_size);
end
end

function q = makeOrthogonalWeights(rows, columns)

[q,r] = qr(randn(rows, columns), 0);
d = diag(r);
q = q * diag(d ./ abs(d));

end
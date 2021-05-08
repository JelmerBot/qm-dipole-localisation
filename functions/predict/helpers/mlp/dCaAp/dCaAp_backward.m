function dLossdX = dCaAp_backward(~, dLossdZ, X)
% dCaAp_backward   Perform backpropagation for a dCaAp layer
%
% Inputs:
% Z - The output from the dCaAp layer.
% dLossdZ - The derivative of the loss function with respect to the output
% of the sigmoid layer. A (H)x(W)x(C)x(N) array.
% X - The input to the dCaAp layer. A (H)x(W)x(C)x(N) array.
%
% Output:
% dLossdX - The derivative of the loss function with respect to the input
% of the dCaAp layer. A (H)x(W)x(C)x(N) array.

%   Copyright 2015-2016 The MathWorks, Inc.

dZdX = dCaAp_derivative(X);
dLossdX = dLossdZ .* dZdX;
end

function y = dCaAp_gpu(x)
    x = gpuarray(x);
    y = 4 .* sigmoid_derivative(x) .* (x > 0);
end

function y = sigmoid_derivative(x)
    t = sigmoid(x);
    y = t .* (1 - t);
end

function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end
function y = dCaAp_derivative(x)
    y = 4 .* sigmoid_second_derivative(x) .* (x > 0);
end

function y = sigmoid_second_derivative(x)
    t = sigmoid(x);
    y = t .* (1 - t) .* (1 - 2 .* t);
end

function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end
function x = relu(x)
% x(x<0) = 0;
x = max(0, x);
end
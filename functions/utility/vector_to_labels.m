function labels = vector_to_labels(vector, prefix, postfix)
    labels = arrayfun(@(x) [prefix, num2str(x), postfix], vector, 'UniformOutput', false);
end
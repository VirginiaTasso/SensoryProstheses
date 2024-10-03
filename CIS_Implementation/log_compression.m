function compressed = log_compression(x, eps)
    eps = 1.0; % to avoid to make the logarithm of 0
    compressed = log10(x + eps);
end
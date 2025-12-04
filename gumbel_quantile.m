function q = gumbel_quantile(alpha)
    % Check that alpha is in the valid range
    assert(alpha > 0 && alpha < 1, 'alpha must be between 0 and 1');

    % Calculate the alpha quantile of the standard Gumbel distribution
    q = -log(-log(alpha));
end
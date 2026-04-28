function tf = isRegularStochastic(P)
    % Returns true if P is column-stochastic AND some power P^k has
    % strictly positive entries (the textbook definition of regular).
    % Searches up to k = n^2 + 1.

    if isa(P, 'sym')
        try
            P = double(P);
        catch
            error('isRegularStochastic: cannot convert symbolic P with free variables to numeric.');
        end
    end

    [m, n] = size(P);
    tol = 1e-9;
    nonneg  = all(P(:) >= -tol);
    sumsOK  = all(abs(sum(P, 1) - 1) <= tol);
    square  = (m == n);

    if ~(nonneg && sumsOK && square)
        fprintf('P is not stochastic, hence not regular.\n');
        if nargout >= 1, tf = false; end
        return
    end

    Pk = P;
    for k = 1:(n^2 + 1)
        if all(Pk(:) > tol)
            fprintf('P is regular: P^%d has all entries > 0.\n', k);
            if nargout >= 1, tf = true; end
            return
        end
        Pk = Pk * P;
    end

    fprintf('P is stochastic but no power up to k = %d has all positive entries (not regular).\n', n^2 + 1);
    if nargout >= 1, tf = false; end
end

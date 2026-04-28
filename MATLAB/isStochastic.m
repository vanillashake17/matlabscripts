function tf = isStochastic(P)
    % Check if P is a (column) stochastic matrix:
    % all entries >= 0 and every column sums to 1.

    if isa(P, 'sym')
        try
            P = double(P);
        catch
            error('isStochastic: cannot convert symbolic P with free variables to numeric.');
        end
    end

    [m, n] = size(P);
    tol = 1e-9;

    nonneg = all(P(:) >= -tol);
    colSums = sum(P, 1);
    sumsOK = all(abs(colSums - 1) <= tol);
    square = (m == n);

    tflocal = nonneg && sumsOK && square;

    if tflocal
        fprintf('P is column-stochastic (size %dx%d).\n', m, n);
    else
        fprintf('P is NOT column-stochastic.\n');
        if ~square,  fprintf('  - not square (size %dx%d)\n', m, n); end
        if ~nonneg,  fprintf('  - has negative entries\n'); end
        if ~sumsOK,  fprintf('  - column sums = %s (should all be 1)\n', mat2str(colSums, 4)); end
    end

    if nargout >= 1, tf = tflocal; end
end

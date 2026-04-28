function B = extendToBasis(V)
    % Given a matrix V whose columns are linearly independent vectors in R^n,
    % extend them to a basis of R^n by appending standard basis vectors.

    [n, k] = size(V);
    if rank(V) < k
        error('extendToBasis: input columns are not linearly independent.');
    end
    if k > n
        error('extendToBasis: more columns than dimension; cannot extend.');
    end

    % Always promote to sym so machine-epsilon doesn't shift pivot selection.
    M = sym([V, eye(n)]);
    pivots = pivot_columns(rref(M));
    Blocal = M(:, pivots);

    addedIdx = pivots(pivots > k) - k;

    fprintf('Extended to a basis of R^%d (added %d standard basis vector(s)).\n', n, n - k);
    if ~isempty(addedIdx)
        fprintf('Appended: e_%s\n', mat2str(addedIdx));
    end
    disp(Blocal);

    if nargout >= 1, B = Blocal; end
end

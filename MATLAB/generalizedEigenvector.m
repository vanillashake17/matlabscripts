function v2 = generalizedEigenvector(A, lambda, v1)
    % Find a generalized eigenvector v2 satisfying (A - lambda I) v2 = v1,
    % where v1 is an eigenvector of A for the (repeated) eigenvalue lambda.
    %
    % Free variables are set to 0 to give one particular v2.

    A = sym(A); v1 = sym(v1(:)); lambda = sym(lambda);
    n = size(A, 1);
    M = A - lambda * eye(n);

    if ~all(simplify(M * v1) == sym(zeros(n, 1)))
        warning('v1 does not satisfy (A - lambda I) v1 = 0; check that v1 is an eigenvector for lambda.');
    end

    aug = [M, v1];
    R = rref(aug);
    pivots = pivot_columns(R);

    if any(pivots == n + 1)
        fprintf('No generalized eigenvector exists: (A - lambda I) v2 = v1 is inconsistent.\n');
        if nargout >= 1, v2 = []; end
        return
    end

    v2local = sym(zeros(n, 1));
    for i = 1:numel(pivots)
        col = pivots(i);
        if col <= n
            v2local(col) = R(i, end);
        end
    end

    fprintf('Generalized eigenvector v2 (free variables set to 0):\n');
    disp(v2local);

    if ~all(simplify(M * v2local - v1) == sym(zeros(n, 1)))
        warning('Verification failed: (A - lambda I) v2 != v1.');
    else
        fprintf('Verified: (A - lambda I) v2 = v1.\n');
    end

    if nargout >= 1, v2 = v2local; end
end


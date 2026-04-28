function W = orthogonalComplement(V)
    % Given a matrix V whose columns span a subspace V of R^n,
    % return a matrix W whose columns are a basis of V^perp = Null(V').
    %
    % Standard fact: V^perp = Null(V^T).

    % Always promote to symbolic so the basis stays exact (no machine-eps
    % leakage from the floating-point null('r') path).
    V = sym(V);
    Wlocal = null(V.');

    [n, ~] = size(V);
    fprintf('V is in R^%d.  dim(V^perp) = %d.\n', n, size(Wlocal, 2));
    fprintf('Basis for V^perp (columns):\n');
    disp(Wlocal);

    if nargout >= 1, W = Wlocal; end
end

function p = charPoly(A)
    % Print the characteristic polynomial det(xI - A), its factorisation,
    % and a table of each eigenvalue with algebraic and geometric multiplicities.

    syms x
    A = sym(A);
    n = size(A, 1);

    plocal = expand(det(x * eye(n) - A));
    fprintf('char(A) = det(xI - A) =\n');
    disp(plocal);

    fprintf('Factored:\n');
    try
        disp(factor(plocal, x));
    catch
        disp(factor(plocal));
    end

    eigs = solve(plocal == 0, x);
    if isempty(eigs)
        fprintf('No symbolic roots found.\n');
        if nargout >= 1, p = plocal; end
        return
    end
    eigs = unique(eigs);

    fprintf('\n  lambda        alg.mult   geo.mult\n');
    fprintf('  ------------  ---------  ---------\n');
    for i = 1:numel(eigs)
        lam = eigs(i);

        % Algebraic multiplicity via repeated derivatives at lambda
        am = 0;
        deriv = plocal;
        while simplify(subs(deriv, x, lam)) == 0
            am = am + 1;
            deriv = diff(deriv, x);
        end

        % Geometric multiplicity = dim Null(A - lam I)
        gm = size(null(A - lam * eye(n)), 2);

        fprintf('  %-12s  %-9d  %-9d\n', char(lam), am, gm);
    end

    if nargout >= 1, p = plocal; end
end

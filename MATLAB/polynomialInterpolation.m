function [coeffs, polyExpr] = polynomialInterpolation(xs, ys, deg)
    % Find a polynomial p(x) = a_deg x^deg + ... + a_1 x + a_0 of degree
    % at most `deg` that satisfies p(xs(i)) = ys(i) for each i.
    %
    % If deg is omitted, the natural choice deg = length(xs) - 1 is used.
    %
    % If the system has a unique solution, returns the numeric/sym coefficients.
    % If it has infinitely many, free coefficients are introduced as t1, t2, ...
    % If it has no solution, errors.
    %
    % Output ordering: coeffs(1) is a_deg (highest), coeffs(end) is a_0.

    if nargin < 3, deg = numel(xs) - 1; end
    xs = sym(xs(:)); ys = sym(ys(:));
    n  = numel(xs);

    % Vandermonde-style matrix: V(i, j) = xs(i)^(deg+1-j)
    V = sym(zeros(n, deg+1));
    for j = 1:deg+1
        V(:, j) = xs.^(deg+1-j);
    end

    Raug = rref([V, ys]);
    pivA = pivot_columns(Raug);
    rA = numel(pivA);
    rV = sum(any(Raug(:, 1:end-1) ~= 0, 2));

    if rA > rV
        error('No polynomial of degree <= %d satisfies these conditions.', deg);
    end

    pivCols  = pivA(pivA <= deg+1);
    freeCols = setdiff(1:deg+1, pivCols);
    coeffs   = sym(zeros(deg+1, 1));

    if isempty(freeCols)
        for i = 1:numel(pivCols)
            coeffs(pivCols(i)) = Raug(i, end);
        end
        fprintf('Unique polynomial.\n');
    else
        params = sym('t', [1, numel(freeCols)]);
        for k = 1:numel(freeCols)
            coeffs(freeCols(k)) = params(k);
        end
        for i = 1:numel(pivCols)
            row = Raug(i, :);
            val = row(end);
            for k = 1:numel(freeCols)
                val = val - row(freeCols(k)) * params(k);
            end
            coeffs(pivCols(i)) = val;
        end
        fprintf('Family of polynomials with %d free parameter(s): %s\n', ...
            numel(freeCols), char(params));
    end

    coeffs = simplify(coeffs);
    syms x;
    polyExpr = simplify(sum(coeffs .* x.^((deg:-1:0)')));

    fprintf('Coefficients [a_%d ; ... ; a_0]:\n', deg); disp(coeffs);
    fprintf('p(x) = %s\n', char(polyExpr));
end

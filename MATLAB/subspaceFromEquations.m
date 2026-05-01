function [basis, dimV] = subspaceFromEquations(C)
% subspaceFromEquations - basis + general solution for V = {x : C*x = 0}.
%
% Given the coefficient matrix C of a homogeneous linear system, returns a
% basis (columns) of the solution subspace and its dimension. Uses exact
% symbolic null space, so output stays in fractions/surds. Also prints the
% general solution to C*x = 0 in the form
%
%     x = s*v1 + t*v2 + r*v3 + t1*v4 + t2*v5 + ...
%
% (free parameters auto-named s, t, r, t1, t2, ... matching least_squares).
%
% Useful for exam questions of the form:
%   "Let V = {(x1, x2, x3, x4) : 5x1 + 3x2 - 2x3 + 3x4 = 0}, find a basis."
% which translates to subspaceFromEquations([5 3 -2 3]).
%
% SYNTAX:
%   [basis, dimV] = subspaceFromEquations(C)
%
% INPUT:
%   C - m-by-n coefficient matrix; each row is one homogeneous equation.
%
% OUTPUT:
%   basis - n-by-dimV matrix whose columns form a basis for V = null(C).
%   dimV  - dim(V) = n - rank(C).

    if ~isnumeric(C) && ~isa(C, 'sym')
        error('Input C must be numeric or symbolic.');
    end
    Csym = sym(C);
    [m, n] = size(Csym);

    fprintf('--- Subspace defined by C*x = 0 ---\n');
    fprintf('Constraint matrix C (size %dx%d):\n', m, n);
    disp(Csym);

    raw_basis = null(Csym);
    dimV = size(raw_basis, 2);

    if dimV == 0
        fprintf('Trivial subspace: only x = 0 satisfies C*x = 0.\n');
        basis = raw_basis;
        return;
    end

    % Raw null() basis: column j has 1 at the j-th free-variable row and 0 at
    % the other free rows. We use this for the literal "x_{free} = param"
    % mapping in the printed general solution. The returned basis is then
    % rescaled (clear_fractions) for cleaner integer display.
    basis = raw_basis;
    for k = 1:dimV
        basis(:, k) = clear_fractions(basis(:, k));
    end

    % Identify free-variable columns from rref(C) so we can label each
    % parameter with the variable it represents.
    R = rref(Csym);
    pivots = pivot_columns(R);
    free_cols = setdiff(1:n, pivots);

    fprintf('rank(C) = %d,   dim(V) = n - rank(C) = %d   (V is in R^%d)\n', ...
        n - dimV, dimV, n);
    fprintf('Basis for V (columns, integer-scaled):\n');
    disp(basis);

    names = param_names(dimV);
    fprintf('Free variables: ');
    for k = 1:dimV
        if k > 1
            fprintf(', ');
        end
        fprintf('x_%d = %s', free_cols(k), names{k});
    end
    fprintf('\n');

    % General solution row-by-row, built from the RAW null basis so that
    % x_{free_cols(k)} = names{k} literally. simplify() collapses the trivial
    % (1*s) entries.
    x_gen = sym(zeros(n, 1));
    for k = 1:dimV
        x_gen = x_gen + sym(names{k}) * raw_basis(:, k);
    end
    x_gen = simplify(x_gen);
    fprintf('General solution to C*x = 0:\n');
    for i = 1:n
        fprintf('   x_%d = %s\n', i, char(x_gen(i)));
    end
    fprintf('   (or in vector form: x = ');
    for j = 1:dimV
        if j > 1
            fprintf(' + ');
        end
        fprintf('%s*%s', names{j}, vec_to_str(raw_basis(:, j)));
    end
    fprintf(')\n');
end

% -------------------------------------------------------------------------
function names = param_names(k)
    base = {'s', 't', 'r'};
    names = cell(1, k);
    for j = 1:k
        if j <= numel(base)
            names{j} = base{j};
        else
            names{j} = sprintf('t%d', j - numel(base));
        end
    end
end

% -------------------------------------------------------------------------
function s = vec_to_str(v)
    parts = strings(numel(v), 1);
    for i = 1:numel(v)
        parts(i) = string(char(v(i)));
    end
    s = char("[" + strjoin(parts, "; ") + "]");
end

% -------------------------------------------------------------------------
function w = clear_fractions(v)
% Scale a symbolic column vector so its entries have no rational denominators
% and share no common integer factor. Basis vectors are defined only up to
% nonzero scalar so this preserves their meaning.
    w = v;
    if isempty(v)
        return;
    end
    try
        if all(isAlways(v == 0))
            return;
        end
    catch
        return;
    end
    try
        [~, D] = numden(v);
        L = D(1);
        for k = 2:numel(D)
            L = lcm(L, D(k));
        end
        scaled = simplify(v * L);

        N = numden(scaled);
        g = N(1);
        for k = 2:numel(N)
            g = gcd(g, N(k));
        end
        if ~isAlways(g == 0) && ~isAlways(g == 1) && ~isAlways(g == -1)
            scaled = simplify(scaled / g);
        end

        w = scaled;
    catch
        w = v;
    end
end

function V = vandermondeMatrix(xs, deg)
% vandermondeMatrix - build a Vandermonde matrix for given x-values up to a
% specified maximum polynomial degree.
%
% For x-values x_1, ..., x_n and max degree d, returns the n-by-(d+1) matrix
%
%      [ x_1^d   x_1^(d-1)   ...   x_1   1 ]
%      [ x_2^d   x_2^(d-1)   ...   x_2   1 ]
%      [   :        :         :     :    : ]
%      [ x_n^d   x_n^(d-1)   ...   x_n   1 ]
%
% i.e. V(i, j) = xs(i) ^ (d + 1 - j). Matches the convention used by
% polynomialInterpolation, so the coefficient vector V*c = y satisfies
% c = [a_d; a_{d-1}; ...; a_1; a_0] for the polynomial p(x) = a_d x^d + ... + a_0.
%
% Symbolic by default — `xs` may be numeric or symbolic, output is sym.
%
% SYNTAX:
%   V = vandermondeMatrix(xs, deg)
%   V = vandermondeMatrix(xs)              % deg defaults to numel(xs) - 1
%
% INPUT:
%   xs  - row or column vector of x-values (numeric or sym).
%   deg - non-negative integer, the highest power allowed in the polynomial.
%         (Default: numel(xs) - 1, the smallest deg that makes V square.)
%
% OUTPUT:
%   V - n-by-(deg+1) symbolic Vandermonde matrix.
%
% EXAMPLES:
%   vandermondeMatrix([0 1 2], 2)        % 3x3 (matches MATLAB's vander([0 1 2]))
%   vandermondeMatrix([0 1 2 3], 5)      % 4x6 (under-determined: deg > #pts)
%   syms a; vandermondeMatrix([1 a a^2], 2)   % parametric inputs OK

    if nargin < 1 || isempty(xs)
        error('vandermondeMatrix: need at least one x-value.');
    end
    xs = sym(xs(:));    % column vector, symbolic
    n  = numel(xs);
    if nargin < 2 || isempty(deg)
        deg = n - 1;
    end
    if ~isnumeric(deg) || ~isscalar(deg) || deg < 0 || deg ~= floor(deg)
        error('vandermondeMatrix: deg must be a non-negative integer.');
    end

    V = sym(zeros(n, deg + 1));
    for j = 1:(deg + 1)
        V(:, j) = xs .^ (deg + 1 - j);
    end

    fprintf('Vandermonde matrix: %d x %d   (n = %d points, max degree d = %d)\n', ...
        n, deg + 1, n, deg);
    fprintf('Convention: V(i, j) = xs(i)^(d + 1 - j)  -> column 1 is x^d, last column is 1.\n');
    disp(V);
end

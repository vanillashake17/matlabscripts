function [inFlag, coeffs] = inSpan(v, S)
    % Decide whether v lies in span(columns of S), using exact symbolic
    % arithmetic. If yes, also return one valid coefficient vector c with
    % S*c = v.
    %
    % Usage:
    %   inSpan(v, S)              -> prints verdict (and coefficients)
    %   tf = inSpan(v, S)         -> just the boolean
    %   [tf, c] = inSpan(v, S)    -> boolean and coefficients

    if ~isa(S, 'sym'), S = sym(S); end
    if ~isa(v, 'sym'), v = sym(v); end
    v = v(:);

    [m, k] = size(S);
    if length(v) ~= m
        error('Dimension mismatch: v has %d entries, S has %d rows.', length(v), m);
    end

    pivS = pivot_columns(rref(S));
    Raug = rref([S, v]);
    pivA = pivot_columns(Raug);
    rS = numel(pivS);
    rA = numel(pivA);

    if rA > rS                 % v introduced a new pivot in last column
        inFlag = false;
        coeffs = [];
        fprintf('v is NOT in span(S).\n');
        return;
    end

    inFlag = true;
    % Particular solution: set free vars = 0, read pivots from last col
    coeffs = sym(zeros(k, 1));
    for idx = 1:rS
        coeffs(pivS(idx)) = Raug(idx, end);
    end

    fprintf('v IS in span(S).\n');
    if nargout < 2
        fprintf('One valid coefficient vector c (so that S*c = v):\n');
        disp(coeffs);
    end
end

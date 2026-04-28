function coords = coordinatesInBasis(v, B)
    % Compute the coordinate vector [v]_B given a basis B (columns) and
    % a vector v in span(B). Uses exact symbolic arithmetic.
    %
    % Errors if columns of B are linearly dependent or if v is not in span(B).

    if ~isa(B, 'sym'), B = sym(B); end
    if ~isa(v, 'sym'), v = sym(v); end
    v = v(:);

    [m, k] = size(B);
    if length(v) ~= m
        error('Dimension mismatch: v has %d entries, B has %d rows.', length(v), m);
    end

    pivB = pivot_columns(rref(B));
    if numel(pivB) < k
        error('Columns of B are linearly dependent -- not a basis.');
    end

    Raug = rref([B, v]);
    pivA = pivot_columns(Raug);
    if numel(pivA) > numel(pivB)
        error('v is not in span(B); cannot represent.');
    end

    coords = sym(zeros(k, 1));
    for i = 1:k
        coords(i) = Raug(i, end);
    end
    coords = simplify(coords);

    fprintf('Coordinates [v]_B (so that B * [v]_B = v):\n');
    disp(coords);
end

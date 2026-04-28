function [RowBasis, ColBasis] = rowColSpace(A)
    % Prints (and optionally returns) bases for Row(A) and Col(A).
    %
    %   rowColSpace(A)              % just displays the bases
    %   [Rb, Cb] = rowColSpace(A)   % also returns them
    %
    % Row space basis = nonzero rows of RREF(A).
    % Col space basis = ORIGINAL columns of A at the pivot positions of RREF(A).
    %
    % Works for both numeric and symbolic matrices.

    % Promote to symbolic so the RREF (and hence the pivot pattern) is exact.
    A = sym(A);
    [m, n] = size(A);
    R = rref(A);
    pivotCols = pivot_columns(R);
    r = numel(pivotCols);

    RowBasisLocal = R(1:r, :);          % each ROW is a basis vector
    ColBasisLocal = A(:, pivotCols);    % each COLUMN is a basis vector

    fprintf('rank(A) = %d   (size %dx%d)\n', r, m, n);
    fprintf('Pivot columns of RREF: %s\n\n', mat2str(pivotCols));

    fprintf('Basis for Row(A)  -- %d vector(s) in R^%d (shown as rows):\n', r, n);
    disp(RowBasisLocal);

    fprintf('Basis for Col(A)  -- %d vector(s) in R^%d (shown as columns):\n', r, m);
    disp(ColBasisLocal);

    if nargout >= 1, RowBasis = RowBasisLocal; end
    if nargout >= 2, ColBasis = ColBasisLocal; end
end

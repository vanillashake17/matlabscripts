function result = isLinearlyIndependent(A)
    % Determine if the columns of A form a linearly independent set,
    % using EXACT symbolic arithmetic (no floating-point tolerance).
    %
    % Input: A - matrix where each column is a vector
    %        (numeric input is auto-converted to sym)
    %
    % Output: true / false (also printed with reason / dependence relation)
    %
    % Edge cases handled:
    %   * empty set            -> independent (vacuously)
    %   * any zero column      -> dependent (reported)
    %   * duplicate columns    -> dependent (caught by rank)
    %   * more vectors than
    %     ambient dimension    -> dependent (caught by rank)
    %   * single non-zero vec  -> independent

    if ~isa(A, 'sym')
        A = sym(A);
    end

    [m, k] = size(A);
    fprintf('Testing %d vector(s) in R^%d.\n', k, m);

    % Empty set: vacuously independent
    if k == 0
        fprintf('Empty set is (vacuously) linearly independent.\n');
        result = true;
        return;
    end

    % Edge case 1: any zero column makes the set dependent
    zeroCols = [];
    for i = 1:k
        if isAlways(all(A(:,i) == 0), 'Unknown', 'false')
            zeroCols(end+1) = i; %#ok<AGROW>
        end
    end
    if ~isempty(zeroCols)
        fprintf('Set contains the zero vector at column(s): %s\n', mat2str(zeroCols));
        fprintf('=> Linearly DEPENDENT (a set containing 0 is always dependent).\n');
        if k >= 2
            other = setdiff(1:k, zeroCols(1));
            fprintf('   e.g. 1*v_%d + 0*(others) = 0  with coefficient on v_%d non-zero.\n', ...
                zeroCols(1), zeroCols(1));
            if isempty(other), end %#ok<NASGU>
        end
        result = false;
        return;
    end

    % Edge case 2: more vectors than ambient dimension
    if k > m
        fprintf('More vectors (%d) than ambient dimension (%d) => DEPENDENT.\n', k, m);
        % Fall through to rank check anyway to print a relation
    end

    % General test: rank via RREF
    R = rref(A);
    pivotCols = pivot_columns(R);
    r = numel(pivotCols);
    fprintf('rank(A) = %d   (need %d for independence)\n', r, k);

    if r == k
        fprintf('=> Linearly INDEPENDENT.\n');
        result = true;
        return;
    end

    % Dependent: build an explicit non-trivial relation using a free column
    freeCols = setdiff(1:k, pivotCols);
    j = freeCols(1);                       % a dependent column
    % R(:, j) gives the coefficients on the pivot columns expressing v_j
    coeffs = sym(zeros(1, k));
    for idx = 1:r
        coeffs(pivotCols(idx)) = R(idx, j);
    end
    coeffs(j) = -1;
    % Relation: sum_i coeffs(i) * v_i = 0
    fprintf('=> Linearly DEPENDENT.\n');
    fprintf('Non-trivial relation (coefficients on v_1..v_%d):\n', k);
    disp(coeffs);
    fprintf('i.e.  ');
    first = true;
    for i = 1:k
        c = coeffs(i);
        if isAlways(c == 0, 'Unknown', 'false'), continue; end
        if first
            fprintf('(%s)*v_%d', char(c), i);
            first = false;
        else
            fprintf(' + (%s)*v_%d', char(c), i);
        end
    end
    fprintf(' = 0\n');
    result = false;
end

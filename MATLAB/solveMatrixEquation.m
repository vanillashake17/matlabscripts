function A = solveMatrixEquation(M, B, varargin)
    % Solve A*M = B for A, using exact symbolic arithmetic.
    %
    % Optional name/value:
    %   'NullVec', N   -- columns of N are known to lie in null(A) (so A*N = 0)
    %                     They are appended to M (with zero columns appended to B).
    %
    % After (optional) augmentation, M must have full row rank for A to be unique.
    % We then pick a maximal set of linearly independent columns of M and solve
    % the resulting square system exactly.

    if ~isa(M, 'sym'), M = sym(M); end
    if ~isa(B, 'sym'), B = sym(B); end

    Nv = sym([]);
    for i = 1:2:numel(varargin)
        switch varargin{i}
            case 'NullVec'
                Nv = sym(varargin{i+1});
            otherwise
                error('Unknown option: %s', varargin{i});
        end
    end

    if ~isempty(Nv)
        if size(Nv, 1) ~= size(M, 1)
            error('NullVec rows (%d) must match M rows (%d).', size(Nv,1), size(M,1));
        end
        M = [Nv, M];
        B = [sym(zeros(size(B, 1), size(Nv, 2))), B];
        fprintf('Augmented with %d null-space column(s).\n', size(Nv, 2));
    end

    m = size(M, 1);
    if rank(M) < m
        error('M does not have full row rank (%d); A is not uniquely determined.', m);
    end

    % pivot_columns derives pivots from the symbolic RREF directly,
    % so we don't lose precision via double(M) on parametric inputs.
    piv = pivot_columns(rref(M));
    pivM = piv(piv <= size(M, 2));
    pivM = pivM(1:m);                 % first m independent columns
    fprintf('Using columns %s of M (a %dx%d invertible block).\n', ...
        mat2str(pivM), m, m);

    Msub = M(:, pivM);
    Bsub = B(:, pivM);

    A = simplify(Bsub / Msub);        % = Bsub * inv(Msub)
    fprintf('Recovered A:\n');
    disp(A);

    % Sanity check: verify A*M = B for ALL columns
    resid = simplify(A * M - B);
    if any(resid(:) ~= 0)
        warning('A*M does not match B on all columns. Check for inconsistency.');
        fprintf('A*M - B =\n'); disp(resid);
    end
end

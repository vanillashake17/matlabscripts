function solveLinearSystem_RREF(A, B)
% SOLVELINEARSYSTEM_RREF Solves A*x = B showing the augmented matrix and RREF.
%   Displays intermediate steps for exam working. Symbolic by default — A
%   and B may contain parameters (e.g. syms a b).

    % --- Step 1: Form Augmented Matrix ---
    A = sym(A);
    B = sym(B(:));
    M = [A, B];
    [~, n] = size(A);

    fprintf('==============================================\n');
    fprintf('             LINEAR SYSTEM SOLVER             \n');
    fprintf('==============================================\n');

    fprintf('\n--- Step 1: Form Augmented Matrix [A | b] ---\n');
    disp(M);   % sym disp shows fractions/surds natively

    % --- Step 2: Compute RREF ---
    R = rref(M);
    pivots = pivot_columns(R);

    fprintf('--- Step 2: Reduced Row Echelon Form (RREF) ---\n');
    disp(R);

    % --- Step 3: Analyze Results ---
    fprintf('--- Step 3: Extract Solution ---\n');

    % Inconsistent: pivot in the augmented (last) column.
    if ismember(n + 1, pivots)
        fprintf('Status: Inconsistent System (No Solution).\n');
        fprintf('Reason: There is a pivot in the last column (0 = nonzero).\n\n');
        return;
    end

    rankA = numel(pivots);
    free_cols = setdiff(1:n, pivots);

    if isempty(free_cols)
        fprintf('Status: Unique Solution.\n');
    else
        fprintf('Status: Infinite Solutions (%d free variable(s)).\n', numel(free_cols));
        fprintf('Free Variables: %s\n', sprintf('x%d ', free_cols));
    end

    % Particular solution (free vars = 0).
    xp = sym(zeros(n, 1));
    xp(pivots) = R(1:rankA, end);
    xp = simplify(xp);

    % Null-space basis: one column per free var.
    N = sym(zeros(n, numel(free_cols)));
    for i = 1:numel(free_cols)
        c = free_cols(i);
        N(c, i) = sym(1);
        N(pivots, i) = -R(1:rankA, c);
    end
    N = simplify(N);

    % --- Print Final Equation ---
    fprintf('\n--- Final General Solution ---\n');

    if isempty(free_cols)
        fprintf('x = %s\n\n', vec_to_str(xp));
        return;
    end

    % parametric form  xp + s*n1 + t*n2 + r*n3 + t1*n4 + ...
    param_names = pick_param_names(numel(free_cols));
    sol_str = vec_to_str(xp);
    for k = 1:numel(free_cols)
        col = N(:, k);
        if all(isAlways(col == 0, 'Unknown', 'false')), continue; end
        sol_str = sprintf('%s + %s*%s', sol_str, param_names{k}, vec_to_str(col));
    end
    fprintf('x = %s\n', sol_str);
    fprintf('   (free parameters %s correspond to free vars x%s)\n\n', ...
        strjoin(param_names, ', '), strjoin(string(free_cols), ', x'));
end

% -------------------------------------------------------------------------
function s = vec_to_str(v)
% Render a column sym vector as "[a; b; c]" with each entry simplified.
    parts = strings(numel(v), 1);
    for i = 1:numel(v)
        parts(i) = string(char(simplify(v(i))));
    end
    s = char("[" + strjoin(parts, "; ") + "]");
end

% -------------------------------------------------------------------------
function names = pick_param_names(k)
% Free parameters auto-named s, t, r, t1, t2, ... (matches least_squares,
% subspaceFromEquations, classify_linear_system).
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

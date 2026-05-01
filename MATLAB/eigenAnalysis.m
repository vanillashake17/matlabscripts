function [P_exact, D_exact] = eigenAnalysis(A, force_orthogonal_if_symmetric)
% eigenAnalysis - Eigenvalue/eigenspace analysis with exact symbolic
% diagonalization, falling back to Jordan form when the matrix is defective.
%
% When every geometric multiplicity equals its algebraic multiplicity, the
% second output D_exact is diagonal and A = P*D*P^{-1}.
%
% When some eigenvalue has gm < am, the function builds Jordan chains by
% solving (A - lambda I) v_{k+1} = v_k and returns the Jordan form J in
% place of D, so that A = P*J*P^{-1}.
%
% SYNTAX:
%   [P_exact, D_exact] = eigenAnalysis(A)
%   [P_exact, D_exact] = eigenAnalysis(A, force_orthogonal_if_symmetric)

% --- 1. Input Validation and Setup ---
if ~isnumeric(A) && ~isa(A, 'sym')
    error('Input must be a numeric or symbolic matrix.');
end
if size(A,1) ~= size(A,2)
    error('Input matrix must be square.');
end
if ~license('test', 'Symbolic_Toolbox')
    error('This function requires the Symbolic Math Toolbox.');
end

% Handle the optional argument. Default to 'true' if not provided.
if nargin < 2
    force_orthogonal_if_symmetric = true;
end

fprintf('--- Starting Symbolic Diagonalization Analysis ---\n');
A_sym = sym(A);
n = size(A_sym, 1);
I = eye(n);
syms x;

% --- 2. Check for Symmetry and Determine Diagonalization Type ---
% [FIX]: Use all(all(...)) to ensure we get a single scalar logical result.
% isAlways returns a matrix of booleans; we need to know if ALL elements match.
check_matrix = isAlways(A_sym == transpose(A_sym));
is_symmetric = all(check_matrix(:)); 

if is_symmetric
    if force_orthogonal_if_symmetric
        fprintf('Matrix is symmetric. Performing ORTHOGONAL diagonalization.\n');
    else
        fprintf('Matrix is symmetric, but orthogonal diagonalization was DISABLED by the user. Performing STANDARD diagonalization.\n');
    end
else
    fprintf('Matrix is not symmetric. Performing STANDARD diagonalization.\n');
end

% --- 3. Calculate Characteristic Polynomial ---
fprintf('\nStep 1: Characteristic Polynomial\n');
char_poly_expr = det(x*I - A_sym);
factored_poly = factor(char_poly_expr);
fprintf('   - Polynomial Form: det(xI - A) = %s\n', char_poly_expr);
fprintf('   - Factored Form: %s\n', string(factored_poly));

% --- 4. Find Eigenvalues and Analyze Eigenspaces ---
fprintf('\nStep 2: Eigenspace Analysis\n');
if is_symmetric && force_orthogonal_if_symmetric
    fprintf('   (Note: A is symmetric, so eigenvectors from different eigenspaces are\n');
    fprintf('    automatically orthogonal. Gram-Schmidt is applied only within an\n');
    fprintf('    eigenspace of dimension > 1.)\n');
end
% MaxDegree=3 lets solve return explicit (Cardano) cubic roots instead of
% the placeholder root(z^k - ..., z, k). Quartics still fall back to root().
eigenvalues_sym = solve(char_poly_expr, x, 'MaxDegree', 3);
try
    eigenvalues_sym = simplify(eigenvalues_sym);
catch
end
unique_eigenvalues = unique(eigenvalues_sym);

% Heuristic: when a symbolic value is unreadably long or contains a
% RootOf-style placeholder, prefer a vpa numeric form for display.
is_ugly = @(s) any(contains(string(s(:)), "root(")) || ...
               any(strlength(string(s(:))) > 80);

is_diagonalizable = true;
P_columns = {};
D_diag_elements = sym([]);
super_diag_flags = []; % flag i = 1 means J(i, i+1) = 1 (Jordan superdiagonal)

for i = 1:length(unique_eigenvalues)
    val = unique_eigenvalues(i);

    % Count algebraic multiplicity
    am = sum(isAlways(eigenvalues_sym == val));

    fprintf('----------------------------------------\n');
    if is_ugly(val)
        fprintf('For Eigenvalue x ≈ %s  (vpa, 6 s.f.)\n', string(vpa(val, 6)));
    else
        fprintf('For Eigenvalue x = %s:\n', string(val));
    end

    eigenspace_matrix = val * I - A_sym;
    basis_vectors = null(eigenspace_matrix);
    gm = size(basis_vectors, 2);

    fprintf('   - Algebraic Multiplicity (AM): %d\n', am);
    fprintf('   - Geometric Multiplicity (GM): %d\n', gm);

    final_basis_for_P = basis_vectors; % Start with the standard basis

    % Perform orthogonalization ONLY if the matrix is symmetric AND the user wants it.
    % Symmetric matrices always have gm = am, so no Jordan branch is needed.
    if is_symmetric && force_orthogonal_if_symmetric
        % Clear rational denominators on each basis vector for cleaner working.
        % Eigenvectors are defined only up to a nonzero scalar, so this is safe.
        for c = 1:gm
            basis_vectors(:,c) = clear_fractions(basis_vectors(:,c));
        end

        fprintf('   - Eigenspace Basis (before Gram-Schmidt):\n');
        for c = 1:gm
            fprintf('     a_%d =\n', c);
            if is_ugly(basis_vectors(:,c))
                disp(vpa(basis_vectors(:,c), 6));
            else
                disp(basis_vectors(:,c));
            end
        end

        if gm == 1
            v = basis_vectors(:,1);
            norm_v = simplify(sqrt(v.' * v));
            unit_v = simplify(v / norm_v);
            fprintf('   - Eigenspace is 1-dimensional; no Gram-Schmidt needed within it.\n');
            fprintf('   - Normalising: ||a_1|| = %s,  w_1 = a_1 / ||a_1|| =\n', string(norm_v));
            if is_ugly(unit_v)
                disp(vpa(unit_v, 6));
            else
                disp(unit_v);
            end
            final_basis_for_P = unit_v;
        else
            fprintf('   - Applying Gram-Schmidt within this eigenspace:\n');
            ortho_basis = sym(zeros(size(basis_vectors)));
            ortho_basis(:,1) = basis_vectors(:,1);
            fprintf('     v_1 = a_1\n');
            for k = 2:gm
                v_k = basis_vectors(:,k);
                formula = sprintf('a_%d', k);
                for j = 1:(k-1)
                    vj  = ortho_basis(:,j);
                    num = simplify(basis_vectors(:,k).' * vj);
                    den = simplify(vj.' * vj);
                    coeff = simplify(num / den);
                    v_k = v_k - coeff * vj;
                    formula = sprintf('%s - (%s)*v_%d', formula, string(coeff), j);
                end
                v_k = simplify(v_k);
                ortho_basis(:,k) = v_k;
                fprintf('     v_%d = %s =\n', k, formula);
                if is_ugly(v_k)
                    disp(vpa(v_k, 6));
                else
                    disp(v_k);
                end
            end

            fprintf('   - Normalising each orthogonal vector:\n');
            unit_basis = sym(zeros(size(ortho_basis)));
            for k = 1:gm
                v = ortho_basis(:,k);
                norm_v = simplify(sqrt(v.' * v));
                unit_basis(:,k) = simplify(v / norm_v);
                fprintf('     ||v_%d|| = %s,  w_%d = v_%d / ||v_%d|| =\n', ...
                    k, string(norm_v), k, k, k);
                if is_ugly(unit_basis(:,k))
                    disp(vpa(unit_basis(:,k), 6));
                else
                    disp(unit_basis(:,k));
                end
            end
            final_basis_for_P = unit_basis;
        end
    else
        % Scale each eigenvector to clear rational denominators.
        % Eigenvectors are defined only up to a nonzero scalar, so this is safe.
        for k = 1:size(basis_vectors, 2)
            basis_vectors(:,k) = clear_fractions(basis_vectors(:,k));
        end
        final_basis_for_P = basis_vectors;

        if is_ugly(basis_vectors)
            fprintf('   - Eigenspace Basis (vpa, 6 s.f.):\n');
            disp(vpa(basis_vectors, 6));
        else
            fprintf('   - Eigenspace Basis:\n');
            disp(basis_vectors);
        end
    end

    if gm < am
        is_diagonalizable = false;
        fprintf('   - Verdict: AM > GM. Defective eigenvalue: building Jordan chain(s).\n');

        [chain_cols, block_sizes] = build_jordan_columns( ...
            A_sym, val, n, am, gm, final_basis_for_P);

        % Display chain(s) and update P / J trackers.
        col_offset = 0;
        for b = 1:numel(block_sizes)
            bsz = block_sizes(b);
            fprintf('   - Jordan block of size %d for lambda = %s:\n', bsz, string(val));
            for k = 1:bsz
                v = chain_cols(:, col_offset + k);
                if k == 1
                    fprintf('     v_%d (eigenvector):\n', k);
                else
                    fprintf('     v_%d (generalized, satisfies (A - lambda I) v_%d = v_%d):\n', ...
                        k, k, k-1);
                end
                if is_ugly(v)
                    disp(vpa(v, 6));
                else
                    disp(v);
                end

                % For generalized eigenvectors, also print the general
                % solution form: particular + s_i * (eigenspace basis).
                % Any vector in null(A - lambda*I) can be added to v_k and
                % still satisfy (A - lambda I) v_k = v_{k-1}.
                if k >= 2
                    num_free = size(final_basis_for_P, 2);
                    free_names = free_param_names(num_free);
                    free_syms = sym(free_names);
                    v_general = v;
                    for j = 1:num_free
                        v_general = v_general + free_syms(j) * final_basis_for_P(:, j);
                    end
                    fprintf('     General form (%s free in R): v_%d =\n', ...
                        strjoin(free_names, ', '), k);
                    if is_ugly(v_general)
                        disp(vpa(v_general, 6));
                    else
                        disp(v_general);
                    end
                end

                P_columns{end+1} = v; %#ok<AGROW>
                D_diag_elements = [D_diag_elements; val]; %#ok<AGROW>
                if k < bsz
                    super_diag_flags = [super_diag_flags; 1]; %#ok<AGROW>
                else
                    super_diag_flags = [super_diag_flags; 0]; %#ok<AGROW>
                end
            end
            col_offset = col_offset + bsz;
        end
    else
        fprintf('   - Verdict: AM = GM. Condition met for this eigenvalue.\n');
        for c = 1:gm
            P_columns{end+1} = final_basis_for_P(:, c); %#ok<AGROW>
            D_diag_elements = [D_diag_elements; val]; %#ok<AGROW>
            super_diag_flags = [super_diag_flags; 0]; %#ok<AGROW>
        end
    end
end

fprintf('----------------------------------------\n');

% --- 5. Final Conclusion and Matrix Construction ---
fprintf('\nStep 3: Final Conclusion and Matrix Construction\n');
if is_diagonalizable
    fprintf('   - The matrix IS diagonalizable.\n');

    % Sort eigenvalues/vectors for cleaner output. Symbolic sort() chokes
    % on complex eigenvalues, so order by (real, imag) of the numeric value.
    try
        numeric_D = double(vpa(D_diag_elements, 16));
        [~, sort_idx] = sortrows([real(numeric_D), imag(numeric_D)]);
    catch
        sort_idx = (1:numel(D_diag_elements))';
    end
    sorted_D_elements = D_diag_elements(sort_idx);
    temp_P = horzcat(P_columns{:});
    P_exact = temp_P(:, sort_idx);
    D_exact = diag(sorted_D_elements);

    fprintf('\n   --- Resulting Matrices ---\n');
    if is_ugly(P_exact) || is_ugly(D_exact)
        fprintf('   P (Eigenvectors as columns) ≈\n');
        disp(vpa(P_exact, 6));
        fprintf('   D (Eigenvalues on diagonal) ≈\n');
        disp(vpa(D_exact, 6));
    else
        fprintf('   P (Eigenvectors as columns) =\n');
        disp(P_exact);
        fprintf('   D (Eigenvalues on diagonal) =\n');
        disp(D_exact);
    end

    if is_symmetric && force_orthogonal_if_symmetric
       fprintf('\n   Note: P is an orthonormal matrix, so inv(P) = transpose(P).\n');
    end
else
    fprintf('   - The matrix is NOT diagonalizable. Returning Jordan form J.\n');

    % Do NOT sort: sorting would break Jordan chain ordering.
    P_exact = horzcat(P_columns{:});
    m = numel(D_diag_elements);
    D_exact = diag(D_diag_elements);
    for k = 1:(m-1)
        if super_diag_flags(k) == 1
            D_exact(k, k+1) = sym(1);
        end
    end

    fprintf('\n   --- Resulting Matrices ---\n');
    if is_ugly(P_exact) || is_ugly(D_exact)
        fprintf('   P (eigenvectors + generalized eigenvectors as columns) ≈\n');
        disp(vpa(P_exact, 6));
        fprintf('   J (Jordan form) ≈\n');
        disp(vpa(D_exact, 6));
    else
        fprintf('   P (eigenvectors + generalized eigenvectors as columns) =\n');
        disp(P_exact);
        fprintf('   J (Jordan form) =\n');
        disp(D_exact);
    end

    % Verify A = P * J * inv(P).
    try
        residual = simplify(A_sym * P_exact - P_exact * D_exact);
        if all(isAlways(residual(:) == 0))
            fprintf('\n   Verified: A * P = P * J  (so A = P * J * inv(P)).\n');
        else
            fprintf('\n   WARNING: A * P != P * J. Manual check recommended.\n');
        end
    catch
        fprintf('\n   (Verification skipped: symbolic check failed.)\n');
    end
end
fprintf('--- Analysis Complete ---\n');
end

% -------------------------------------------------------------------------
function [chain_cols, block_sizes] = build_jordan_columns(A_sym, lambda, n, am, gm, eig_basis)
% Build the columns of P for one defective eigenvalue.
%
% Returns chain_cols (n by am) with chains laid out left-to-right:
%   [v1_(1), v1_(2), ..., v1_(k1), v2_(1), v2_(2), ..., v2_(k2), ...]
% where each v_(j+1) satisfies (A - lambda I) v_(j+1) = v_(j).
% block_sizes is a row vector of chain lengths summing to am.

    M = A_sym - lambda * eye(n);

    if gm == 1
        % Common exam case: a single Jordan block of size am.
        % Do NOT clear_fractions on chain elements: rescaling v_k alone
        % breaks the Jordan recurrence (A - lambda I) v_k = v_{k-1}.
        % Leave v_k with whatever fractions solve_chain_step produces so
        % the pairing with v_{k-1} is exact.
        chain = sym(zeros(n, am));
        chain(:, 1) = eig_basis(:, 1);
        for k = 2:am
            chain(:, k) = solve_chain_step(M, chain(:, k-1), n, lambda, k);
        end
        chain_cols = chain;
        block_sizes = am;
        return;
    end

    % General case (1 < gm < am): use null spaces of M^k to determine the
    % block structure, then build chains top-down.
    null_dims = zeros(1, am);
    null_bases = cell(1, am);
    for k = 1:am
        Nk = null(M^k);
        null_bases{k} = Nk;
        null_dims(k) = size(Nk, 2);
    end

    % Number of blocks of size exactly k:  d_k - d_{k-1}  -  (d_{k+1} - d_k)
    % with d_0 = 0 and d_{am+1} = d_am.
    d_seq = [0, null_dims, null_dims(am)];
    block_count = zeros(1, am);
    for k = 1:am
        block_count(k) = (d_seq(k+1) - d_seq(k)) - (d_seq(k+2) - d_seq(k+1));
    end

    chain_cols = sym(zeros(n, 0));
    block_sizes = [];
    used = sym(zeros(n, 0));

    for k = am:-1:1
        for b = 1:block_count(k)
            if k == 1
                exclude = used;
            else
                exclude = [used, null_bases{k-1}];
            end
            v_top = pick_independent_column(null_bases{k}, exclude);

            chain = sym(zeros(n, k));
            chain(:, k) = v_top;
            for j = (k-1):-1:1
                chain(:, j) = M * chain(:, j+1);
            end
            % Do NOT clear_fractions per column here: each chain step
            % satisfies (A - lambda I) v_{j+1} = v_j exactly, and rescaling
            % a single v_j by an LCM would break the recurrence with its
            % neighbour. Leave whatever fractions arise from null(M^k).
            chain_cols = [chain_cols, chain]; %#ok<AGROW>
            block_sizes = [block_sizes, k]; %#ok<AGROW>
            used = [used, chain]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function v = solve_chain_step(M, v_prev, n, lambda, k)
% Solve M * v = v_prev for the next vector in a Jordan chain.
% Free variables are pinned to 0 (any particular solution will do).
    aug = [M, v_prev];
    R = rref(aug);
    pivots = pivot_columns(R);
    if any(pivots == n + 1)
        error('eigenAnalysis:chainBreak', ...
            ['Cannot extend Jordan chain at step %d for lambda = %s: ', ...
             '(A - lambda I) v_%d = v_%d is inconsistent.'], ...
            k, char(string(lambda)), k, k-1);
    end
    v = sym(zeros(n, 1));
    for i = 1:numel(pivots)
        col = pivots(i);
        if col <= n
            v(col) = R(i, end);
        end
    end
end

% -------------------------------------------------------------------------
function w = pick_independent_column(pool, used)
% Return the first column of pool that is linearly independent of used.
    if isempty(used) || size(used, 2) == 0
        w = pool(:, 1);
        return;
    end
    base_rank = double(rank(used));
    for i = 1:size(pool, 2)
        if double(rank([used, pool(:, i)])) > base_rank
            w = pool(:, i);
            return;
        end
    end
    error('eigenAnalysis:noIndepVec', ...
        'Could not find a vector in null(M^k) independent of the existing chain basis.');
end

% -------------------------------------------------------------------------
function names = free_param_names(k)
% Return k free-parameter names following the project convention:
% s, t, r, then t1, t2, t3, ...  (matches least_squares / subspaceFromEquations).
    base = {'s', 't', 'r'};
    names = cell(1, k);
    for i = 1:k
        if i <= numel(base)
            names{i} = base{i};
        else
            names{i} = sprintf('t%d', i - numel(base));
        end
    end
end

% -------------------------------------------------------------------------
function w = clear_fractions(v)
% Scale a symbolic column vector so its entries have no rational
% denominators and share no common integer factor. Eigenvectors are defined
% only up to a nonzero scalar, so this preserves their meaning. Returns v
% unchanged if anything fails (irrational/complex denominators that don't
% admit a rational lcm/gcd, etc.).
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

        % Pull out a common integer factor from the (now denominator-free) entries.
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
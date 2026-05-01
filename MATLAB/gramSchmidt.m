function [U, Q, dependent_indices] = gramSchmidt(V, tolerance)
% gramSchmidt Performs Gram-Schmidt orthogonalisation with full step-by-step
% workings printed as fractions/surds (symbolic internally), and linear
% dependence detection.
%
% Inputs:
%   V - Matrix (numeric or sym) whose columns are the input vectors
%   tolerance - (optional) Numeric threshold for dependence detection
%               (default 1e-10). Used as a fallback after the exact
%               symbolic isAlways(...==0) check fails.
%
% Outputs:
%   U - Orthogonal (unnormalised) vectors, returned as a numeric matrix
%   Q - Orthonormal vectors, returned as a numeric matrix
%   dependent_indices - Indices of linearly dependent columns in V
%
% Saved automatically in workspace as:
%   orthogonal_set_gramSchmidt
%   orthonormal_set_gramSchmidt        (factored struct, exam-style)
%   orthonormal_matrix_gramSchmidt     (raw numeric matrix)
%   dependent_indices_gramSchmidt
%
% Example:
%   V = [1 0 1; 1 1 0; 0 1 1];
%   gramSchmidt(V);

    if nargin < 2
        tolerance = 1e-10;
    end

    % Promote to symbolic internally so workings stay in fraction/surd form.
    % Inputs may be numeric or sym; we always print via char(sym(...)).
    V_sym = sym(V);
    [~, n] = size(V_sym);
    U_sym = sym([]);
    Q_sym = sym([]);
    dependent_indices = [];

    fprintf('\n=== Gram-Schmidt Process ===\n\n');
    for j = 1:n
        fprintf('--- Processing v%d ---\n', j);
        fprintf('  v%d = ', j); print_sym_row(V_sym(:, j).');

        u_j = V_sym(:, j);
        kept = size(U_sym, 2);

        % Subtract projections onto previous orthogonal vectors
        for i = 1:kept
            wu_i  = simplify(dot(V_sym(:, j), U_sym(:, i)));
            uu_i  = simplify(dot(U_sym(:, i), U_sym(:, i)));
            coeff = simplify(wu_i / uu_i);
            proj  = simplify(coeff * U_sym(:, i));

            fprintf('  Projection onto u%d:\n', i);
            fprintf('    <v%d, u%d> = %s\n', j, i, char(wu_i));
            fprintf('    <u%d, u%d> = %s\n', i, i, char(uu_i));
            fprintf('    c%d = <v%d,u%d>/<u%d,u%d> = %s\n', ...
                i, j, i, i, i, char(coeff));
            fprintf('    proj_{u%d}(v%d) = c%d * u%d = ', i, j, i, i);
            print_sym_row(proj.');

            u_j = simplify(u_j - proj);
        end

        if kept > 0
            fprintf('  u%d = v%d - (sum of projections) = ', j, j);
            print_sym_row(u_j.');
        else
            fprintf('  u%d = v%d (no previous u_i to subtract)\n', j, j);
        end

        nrm_uj_sym = simplify(sqrt(dot(u_j, u_j)));
        fprintf('  ||u%d|| = %s\n', j, char(nrm_uj_sym));

        % Dependence check: exact symbolic zero, else numeric tolerance fallback
        if isAlways(nrm_uj_sym == 0)
            is_dependent = true;
        else
            try
                is_dependent = double(nrm_uj_sym) < tolerance;
            catch
                is_dependent = false;
            end
        end

        if is_dependent
            fprintf('  -> ||u%d|| < tolerance (%.0e). v%d is LINEARLY DEPENDENT (ignored).\n\n', ...
                j, tolerance, j);
            dependent_indices = [dependent_indices, j];
            continue;
        end

        U_sym = [U_sym, u_j];
        q_j = simplify(u_j / nrm_uj_sym);
        Q_sym = [Q_sym, q_j];

        fprintf('  q%d = u%d / ||u%d|| = ', j, j, j);
        print_sym_row(q_j.');
        fprintf('  -> v%d added to orthogonal and orthonormal sets.\n\n', j);
    end

    % Convert sym matrices back to double for legacy returns / factored display
    U = double(U_sym);
    Q = double(Q_sym);

    % Display results
    fprintf('=== Summary ===\n');
    fprintf('Independent vectors: %d\n', size(Q, 2));
    if isempty(dependent_indices)
        fprintf('All vectors are linearly independent!\n\n');
    else
        fprintf('Dependent vector indices: [ ');
        fprintf('%d ', dependent_indices);
        fprintf(']\n\n');
    end

    fprintf('=== Orthogonal Set (U) ===\n');
    disp(U_sym);
    fprintf('=== Orthonormal Set (Q) - Factored Form ===\n');
    Q_factored = build_factored_struct(U);
    print_factored_struct(Q_factored);
    fprintf('\n=== Verification (Q^T * Q) ===\n');
    disp(simplify(Q_sym.' * Q_sym));

    % Assign results to base workspace
    assignin('base', 'orthogonal_set_gramSchmidt', U);
    assignin('base', 'orthonormal_set_gramSchmidt', Q_factored);  % struct with scalar + int_vec
    assignin('base', 'orthonormal_matrix_gramSchmidt', Q);        % raw numeric matrix (for computation)
    assignin('base', 'dependent_indices_gramSchmidt', dependent_indices);

    fprintf('\nResults saved to workspace:\n');
    fprintf('  orthogonal_set_gramSchmidt       - orthogonal vectors (matrix)\n');
    fprintf('  orthonormal_set_gramSchmidt      - factored form struct (scalar + integer vector)\n');
    fprintf('  orthonormal_matrix_gramSchmidt   - orthonormal vectors as numeric matrix\n');
    fprintf('  dependent_indices_gramSchmidt\n');
    fprintf('\nTo display orthonormal_set_gramSchmidt:\n');
    fprintf('  print_factored_struct(orthonormal_set_gramSchmidt)\n');
    fprintf('=============================================\n\n');
end


function print_sym_row(row)
% Print a symbolic row vector in [a, b, c] form (preserves fractions/surds).
    fprintf('[');
    for k = 1:length(row)
        if k == 1
            fprintf('%s', char(row(k)));
        else
            fprintf(', %s', char(row(k)));
        end
    end
    fprintf(']\n');
end


function Q_factored = build_factored_struct(U)
% Builds a struct array representing each orthonormal vector in factored form:
%   q_k = (1 / scalar_val) * int_vec
%
% Each struct entry has fields:
%   .scalar_val  - numeric scalar (the denominator)
%   .scalar_str  - string representation of scalar (e.g. '2*sqrt(11)')
%   .int_vec     - integer column vector (numerator)
%   .index       - column index k

    n = size(U, 2);
    Q_factored = struct('index', {}, 'scalar_val', {}, 'scalar_str', {}, 'int_vec', {});

    for k = 1:n
        u = U(:, k);
        nrm_sq = round(dot(u, u) * 1e9) / 1e9;
        nrm = sqrt(nrm_sq);

        [int_vec, denom] = clear_fractions(u);
        scalar_val = denom * nrm;
        scalar_str = format_scalar(scalar_val);

        Q_factored(k).index      = k;
        Q_factored(k).scalar_val = scalar_val;
        Q_factored(k).scalar_str = scalar_str;
        Q_factored(k).int_vec    = int_vec;
    end
end


function print_factored_struct(Q_factored)
% Pretty-prints a Q_factored struct (as returned by build_factored_struct or
% stored in orthonormal_set_gramSchmidt).

    fprintf('\n  { ');
    n = length(Q_factored);
    for k = 1:n
        int_vec    = Q_factored(k).int_vec;
        scalar_str = Q_factored(k).scalar_str;

        fprintf('(1/%s) * [', scalar_str);
        for r = 1:length(int_vec)
            if r == 1
                fprintf('%g', int_vec(r));
            else
                fprintf(', %g', int_vec(r));
            end
        end
        fprintf(']');

        if k < n
            fprintf(',\n    ');
        end
    end
    fprintf(' }\n');
end


function [int_vec, denom] = clear_fractions(v)
% Convert a rational vector to integer form by finding common denominator.
% Returns int_vec (integers) and denom such that v = int_vec / denom.
% The int_vec is reduced by its GCD for simplest form.

    tol = 1e-9;
    % Find denominator for each entry by checking multiples up to max_d
    max_d = 10000;
    denom = 1;
    for i = 1:length(v)
        for d = 1:max_d
            if abs(v(i) * d - round(v(i) * d)) < tol
                denom = lcm(denom, d);
                break;
            end
        end
    end
    int_vec = round(v * denom);

    % Divide out GCD of all entries to get simplest integer form
    g = 0;
    for i = 1:length(int_vec)
        g = gcd(g, abs(int_vec(i)));
    end
    if g > 1
        int_vec = int_vec / g;
        denom = denom / g;
    end
end


function str = format_scalar(val)
% Format a scalar value as a clean string.
% Tries to express val as a/sqrt(b) or sqrt(b)/a form.
%
%   val = p/q * sqrt(r)  =>  displayed as "q/(p*sqrt(r))" etc.

    tol = 1e-6;

    % Check if val is a simple integer or simple fraction
    for q = 1:200
        p = val * q;
        if abs(p - round(p)) < tol
            p = round(p);
            g = gcd(abs(p), q);
            p = p / g;
            q_red = q / g;
            if q_red == 1
                str = sprintf('%d', p);
            else
                str = sprintf('%d/%d', p, q_red);
            end
            return;
        end
    end

    % Check if val = sqrt(b) for integer b
    b = val^2;
    b_r = round(b);
    if abs(b - b_r) < tol
        if b_r == 1
            str = '1';
        else
            str = sprintf('sqrt(%d)', b_r);
        end
        return;
    end

    % Check if val = (p/q) * sqrt(r) for small integers
    for r = 2:500
        test = val / sqrt(r);   % should be p/q
        for q = 1:100
            p = test * q;
            if abs(p - round(p)) < tol && round(p) > 0
                p = round(p);
                g = gcd(p, q);
                p = p/g; q = q/g;
                if q == 1
                    if p == 1
                        str = sprintf('sqrt(%d)', r);
                    else
                        str = sprintf('%d*sqrt(%d)', p, r);
                    end
                else
                    if p == 1
                        str = sprintf('sqrt(%d)/%d', r, q);
                    else
                        str = sprintf('%d*sqrt(%d)/%d', p, r, q);
                    end
                end
                return;
            end
        end
    end

    % Fallback: decimal
    str = sprintf('%.6g', val);
end

function [U, Q, dependent_indices] = gramSchmidt(V, tolerance)
% gramSchmidt Performs Gram-Schmidt orthogonalization with linear dependence detection
%
% Inputs:
%   V - Matrix where each column is a vector
%   tolerance - (optional) Threshold for dependence detection (default 1e-10)
%
% Outputs:
%   U - Orthogonal (unnormalized) vectors
%   Q - Orthonormal vectors
%   dependent_indices - Indices of linearly dependent columns in V
%
% Saved automatically in workspace as:
%   orthogonal_set_gramSchmidt
%   orthonormal_set_gramSchmidt
%   dependent_indices_gramSchmidt
%
% Example:
%   V = [1 0 1; 1 1 0; 0 1 1];
%   gramSchmidt(V);

    if isa(V, 'sym')
        error('gramSchmidt:symInput', ...
            ['gramSchmidt is numeric-only (uses norm() < tolerance). ', ...
             'For symbolic input use orthonormalize_rational(V) ', ...
             '(exact rational/surd output) or orthogonalDiagonalize(A) for symmetric A.']);
    end
    if nargin < 2
        tolerance = 1e-10;
    end

    [~, n] = size(V);
    U = [];
    Q = [];
    dependent_indices = [];

    fprintf('\n=== Gram-Schmidt Process ===\n\n');
    for j = 1:n
        fprintf('Processing vector v%d:\n', j);
        u_j = V(:, j);

        % Subtract projections onto previous orthogonal vectors
        for i = 1:size(U, 2)
            proj = (dot(V(:, j), U(:, i)) / dot(U(:, i), U(:, i))) * U(:, i);
            u_j = u_j - proj;
        end

        % Check for linear dependence
        if norm(u_j) < tolerance
            fprintf('  -> v%d is LINEARLY DEPENDENT (ignored)\n\n', j);
            dependent_indices = [dependent_indices, j];
            continue;
        end

        % Add orthogonal and orthonormal vectors
        U = [U, u_j];
        q_j = u_j / norm(u_j);
        Q = [Q, q_j];

        fprintf('  -> v%d added to orthogonal and orthonormal sets (||u%d|| = %.4f)\n\n', j, j, norm(u_j));
    end

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
    disp(U);
    fprintf('=== Orthonormal Set (Q) - Factored Form ===\n');
    Q_factored = build_factored_struct(U);
    print_factored_struct(Q_factored);
    fprintf('\n=== Verification (Q^T * Q) ===\n');
    disp(Q' * Q);

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

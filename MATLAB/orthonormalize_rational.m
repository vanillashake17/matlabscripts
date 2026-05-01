function E = orthonormalize_rational(U)
% ORTHONORMALIZE_RATIONAL Orthonormalises a set of vectors using exact
% symbolic (rational/surd) Gram-Schmidt, with full workings printed.
%
% Usage:
%   E = orthonormalize_rational(U)
%
% Input:
%   U - matrix whose columns are the vectors v1, v2, ... to orthonormalise.
%
% Output:
%   E - symbolic matrix whose columns are orthonormal vectors w1, w2, ...
%
% Pipeline (printed step-by-step):
%   1) v_k                                     (input vector)
%   2) for each previous w_j:                  print <v_k, w_j>
%      v_k <- v_k - <v_k, w_j> * w_j
%   3) print v_k after subtraction             (orthogonal vector u_k)
%   4) ||u_k||                                  (exact symbolic norm)
%   5) w_k = u_k / ||u_k||                     (orthonormal vector)
%
% NOTE: The original MA1508E convention writes Gram-Schmidt using the
% *unnormalised* orthogonal vectors u_j (not the orthonormal w_j). With
% w_j's the formula reduces to v_k - sum_j <v_k, w_j> w_j because
% <w_j, w_j> = 1, which is exactly what is shown below.
%
% Errors only if a column is linearly dependent on earlier columns
% (i.e. the orthogonal residual u_k collapses to zero).

    U_sym = sym(U);
    [n, k] = size(U_sym);
    E = sym(zeros(n, k));

    fprintf('\n=== Orthonormalisation (exact symbolic Gram-Schmidt) ===\n\n');

    for i = 1:k
        fprintf('--- Processing v%d ---\n', i);
        v0 = U_sym(:, i);
        fprintf('  v%d = ', i); disp(v0.');

        v = v0;
        for j = 1:i-1
            cij = simplify(dot(v0, E(:, j)));   % <v_i, w_j>; w_j is unit
            fprintf('  <v%d, w%d> = %s\n', i, j, char(cij));
            fprintf('  proj_{w%d}(v%d) = <v%d,w%d> * w%d = ', j, i, i, j, j);
            disp((simplify(cij * E(:, j))).');
            v = simplify(v - cij * E(:, j));
        end

        if i > 1
            fprintf('  u%d = v%d - sum of projections = ', i, i);
            disp(v.');
        else
            fprintf('  u%d = v%d (no previous w_j)\n', i, i);
        end

        if isAlways(v == sym(zeros(n,1)))
            error('orthonormalize_rational:dependent', ...
                'Vector v%d is linearly dependent on previous vectors.', i);
        end

        nrm = simplify(sqrt(dot(v, v)));
        fprintf('  ||u%d|| = sqrt(<u%d,u%d>) = %s\n', i, i, i, char(nrm));

        w = simplify(v / nrm);
        fprintf('  w%d = u%d / ||u%d|| = ', i, i, i);
        disp(w.');
        fprintf('\n');

        E(:, i) = w;
    end

    fprintf('=== Orthonormal set (columns of E) ===\n');
    disp(E);

    fprintf('=== Verification (E^T * E should be I) ===\n');
    disp(simplify(E.' * E));
    fprintf('========================================================\n\n');
end

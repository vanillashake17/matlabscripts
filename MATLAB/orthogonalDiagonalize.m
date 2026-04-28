function [P, D] = orthogonalDiagonalize(A)
    % For a symmetric matrix A, return orthogonal P and diagonal D with
    % A = P D P^T (equivalently D = P^T A P).
    % Uses: eigenspace bases -> Gram-Schmidt within each eigenspace -> normalize.

    A = sym(A);
    n = size(A, 1);

    if ~isequal(A, A.')
        error('orthogonalDiagonalize: A must be symmetric.');
    end

    eigs = unique(eig(A));

    Plocal = sym(zeros(n, 0));
    Dlocal = sym([]);

    for i = 1:numel(eigs)
        lam = eigs(i);
        N = null(A - lam * eye(n));         % basis for eigenspace
        Q = local_gramSchmidt(N);            % orthogonalize
        for j = 1:size(Q, 2)                 % normalize
            Q(:, j) = Q(:, j) / sqrt(simplify(Q(:, j).' * Q(:, j)));
        end
        Q = simplify(Q);
        Plocal = [Plocal, Q]; %#ok<AGROW>
        Dlocal = blkdiag(Dlocal, lam * eye(size(Q, 2)));
    end

    Plocal = simplify(Plocal);
    Dlocal = simplify(Dlocal);

    if ~isequal(simplify(Plocal * Dlocal * Plocal.' - A), sym(zeros(n)))
        warning('Verification A = P*D*P^T failed; inspect output.');
    else
        fprintf('Verified: A = P * D * P^T.\n');
    end

    fprintf('Orthogonal P (columns = orthonormal eigenvectors):\n');
    disp(Plocal);
    fprintf('Diagonal D (eigenvalues in matching order):\n');
    disp(Dlocal);

    if nargout >= 1, P = Plocal; end
    if nargout >= 2, D = Dlocal; end
end

function Q = local_gramSchmidt(N)
    [n, k] = size(N);
    Q = sym(zeros(n, 0));
    for j = 1:k
        v = N(:, j);
        for i = 1:size(Q, 2)
            v = v - (Q(:, i).' * v) / (Q(:, i).' * Q(:, i)) * Q(:, i);
        end
        v = simplify(v);
        if any(v ~= 0)
            Q = [Q, v]; %#ok<AGROW>
        end
    end
end

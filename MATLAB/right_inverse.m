function R = right_inverse(A)
%right_inverse Right inverse R of A so that A*R = I (full row rank required).
    % Formula: R = A' * inv(AA')
    % Condition: Rows must be linearly independent (Full Row Rank)
    
    % Always promote to symbolic so machine-epsilon noise doesn't leak in
    % from floating-point inv((AA')).
    A = sym(A);
    [m, ~] = size(A);
    r = rank(A);

    if r == m
        R = simplify(A.' * inv(A * A.'));
        fprintf('[Success] Right Inverse found.\n');
    else
        R = sym([]);
        fprintf('[Flag] NO Right Inverse. Rows are dependent (Rank %d < Rows %d).\n', r, m);
    end
end
function L = left_inverse(A)
%left_inverse Left inverse L of A so that L*A = I (full column rank required).
    % Formula: L = inv(A'A) * A'
    % Condition: Columns must be linearly independent (Full Column Rank)
    
    % Always promote to symbolic so machine-epsilon noise doesn't leak in
    % from floating-point inv((A'A)).
    A = sym(A);
    [~, n] = size(A);
    r = rank(A);

    if r == n
        L = simplify(inv(A.' * A) * A.');
        fprintf('[Success] Left Inverse found.\n');
    else
        L = sym([]);
        fprintf('[Flag] NO Left Inverse. Columns are dependent (Rank %d < Cols %d).\n', r, n);
    end
end
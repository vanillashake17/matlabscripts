function [U, S, V] = svd_exact(A)
% svd_exact - Performs an exact, symbolic Singular Value Decomposition.
%
% [FINAL CORRECTED VERSION] - Replaces the unreliable null() method for
% completing the U matrix with a robust QR decomposition on an augmented
% matrix. This definitively fixes the dimension mismatch errors.
%
% This function requires the Symbolic Math Toolbox and 'eigenAnalysis.m'.
%
% SYNTAX: [U, S, V] = svd_exact(A)

% --- 1. Input Validation and Setup ---
if ~isnumeric(A) && ~isa(A, 'sym')
    error('Input must be a numeric or symbolic matrix.');
end
if ~exist('eigenAnalysis.m', 'file')
    error('This function requires the eigenAnalysis.m function to be in the same path.');
end
if ~license('test', 'Symbolic_Toolbox')
    error('This function requires the Symbolic Math Toolbox.');
end

fprintf('--- Starting Exact Singular Value Decomposition ---\n');
A_sym = sym(A);
[m, n] = size(A_sym);
disp('Input Matrix A:');
disp(A_sym);

% --- 2. Analyze B = A'*A to find V and Singular Values ---
fprintf('\nStep 1: Orthogonally diagonalize B = A''*A\n');
B = A_sym' * A_sym;
[V, D_squared] = eigenAnalysis(B);

eigenvalues = diag(D_squared);
[sorted_eigenvalues, sort_idx] = sort(eigenvalues, 'descend');
V = V(:, sort_idx);

fprintf('--- Eigenvalues of A''*A (λ) in descending order ---\n');
disp(sorted_eigenvalues.');
singular_values = simplify(sqrt(sorted_eigenvalues));
fprintf('--- Corresponding Singular Values (σ = sqrt(λ)) ---\n');
disp(singular_values.');
fprintf('--- Resulting Orthonormal Matrix V ---\n');
disp(V);

% --- 3. Construct the Sigma Matrix (S) ---
fprintf('\nStep 2: Construct the Sigma (S) matrix\n');
S = sym(zeros(m, n));
rank_A = sum(isAlways(singular_values > 0));
if rank_A > 0
    S(1:rank_A, 1:rank_A) = diag(singular_values(1:rank_A));
end
disp('Matrix S:');
disp(S);

% --- 4. Construct the U Matrix ---
fprintf('\nStep 3: Construct the Orthonormal Matrix U\n');
U_r = sym(zeros(m, rank_A));
for i = 1:rank_A
    U_r(:, i) = (1/singular_values(i)) * A_sym * V(:, i);
end
U_r = simplify(U_r);

fprintf('--- First r=%d columns of U (calculated via u_i = (1/σ_i)Av_i) ---\n', rank_A);
disp(U_r);

% [FIX] Use the robust augmented matrix + QR method to complete the basis.
if m > rank_A
    fprintf('--- Completing U by extending the basis with QR decomposition ---\n');
    % Augment U_r with the identity matrix to guarantee the columns span R^m.
    temp_basis = [U_r, eye(m)];
    
    % The qr function will create a full orthonormal basis.
    [Q, ~] = qr(temp_basis);
    
    % The final U is the first m columns of the resulting Q.
    U = Q(:, 1:m);
    fprintf('--- Orthonormal basis completed. Final U is %d x %d. ---\n', m, m);
else
    U = U_r;
end
U = simplify(U);
fprintf('--- Final Orthonormal Matrix U ---\n');
disp(U);

% --- 5. Verification ---
fprintf('\nStep 4: Verification Check\n');
A_reconstructed = simplify(U * S * V');
disp('Verifying A = U*S*V'':');
disp(A_reconstructed);

if isAlways(A_reconstructed == A_sym)
    fprintf('\nVerification successful! The reconstructed matrix matches the original A.\n');
else
    fprintf('\nVerification FAILED! The reconstructed matrix does NOT match A.\n');
end
fprintf('--- SVD Complete ---\n');
end
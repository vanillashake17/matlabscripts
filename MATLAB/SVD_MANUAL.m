clear; clc;
fprintf('--- SVD Construction Script ---\n\n');

%% PART 0: Define Matrix A
A = sym([1 0 1 0; -1 2 1 2; 1 1 2 1]);
fprintf('Matrix A:\n'); disp(A);

%% PART 1: The "How we got here" (Automatic Analysis)
% This section calculates the raw eigenvalues and eigenvectors of A'*A
% so you know what numbers to put in the Manual section.

fprintf('--- Part 1: Analyzing A''*A to find Sigmas and v''s ---\n');
ATA = A' * A;
[V_raw, D_raw] = eig(ATA);

% Extract eigenvalues from the diagonal matrix D
eigenvalues = diag(D_raw);

fprintf('Calculated Eigenvalues of A''*A (Squared Sigmas):\n');
disp(eigenvalues);
fprintf('Calculated Eigenvectors (Columns match the values above):\n');
disp(V_raw);
fprintf('NOTE: MATLAB usually sorts Low -> High. You want the High ones.\n');
fprintf('----------------------------------------------------------\n\n');

%% PART 2: Manual Overrides (Define your Sigmas & v vectors here)
% Based on Part 1, pick the non-zero eigenvalues (e.g., 14 and 5).
% Sqrt(14) -> sigma1, Sqrt(5) -> sigma2.
% Pick the corresponding columns from V_raw for v1 and v2.

fprintf('--- Part 2: Defining Parameters ---\n');

% --- Pair 1 (Largest Eigenvalue: 14) ---
sigma1 = sym(sqrt(14)); 
% Normalizing just to be safe, though eig() usually does it.
v1 = sym((1/sqrt(3)) * [0; 1; 1; 1]); 

% --- Pair 2 (Next Eigenvalue: 5) ---
sigma2 = sym(sqrt(5));
v2 = sym((1/sqrt(15)) * [-3; 1; -2; 1]);

fprintf('Using Sigma1: %s\n', char(sigma1));
fprintf('Using Sigma2: %s\n', char(sigma2));

%% PART 3: Compute V (Completes the basis using QR)
fprintf('\n--- Part 3: Computing full V (Right Singular Vectors) ---\n');

% QR automatically finds orthogonal vectors (v3, v4) to fill the space
% We feed it [v1 v2], it gives us [v1 v2 v3 v4] perfectly orthonormal.
[V, ~] = qr([v1, v2]); 
V = simplify(V);

% Extract specific vectors for display if needed
v3 = V(:,3); 
v4 = V(:,4);

fprintf('Full V Matrix:\n'); disp(V);
fprintf('v3 (Null space basis):\n'); disp(v3);
fprintf('v4 (Null space basis):\n'); disp(v4);

%% PART 4: Compute U (Calculate u1, u2 -> Complete with QR)
fprintf('\n--- Part 4: Computing full U (Left Singular Vectors) ---\n');

% Calculate active u vectors using formula: u = (1/sigma) * A * v
u1 = simplify((1/sigma1) * A * v1);
u2 = simplify((1/sigma2) * A * v2);

% QR automatically finds u3 to fill the space for the target R3 dimensions
[U, ~] = qr([u1, u2]);
U = simplify(U);

fprintf('Full U Matrix:\n'); disp(U);

%% PART 5: Verification
fprintf('\n--- Part 5: Verification ---\n');
% Reconstruct Sigma Matrix S
S = sym(zeros(size(A)));
S(1,1) = sigma1; 
S(2,2) = sigma2;

% Check A - USV'
Check = simplify(A - (U * S * V'));
if isequal(Check, sym(zeros(size(A))))
    fprintf('SUCCESS: A matches U * S * V'' exactly.\n');
else
    fprintf('WARNING: Reconstruction error.\n');
    disp(Check);
end

x=[-1;-0.15;0;1;2.15;3;3.3];
y=[0;3;3;0;-3;0;3];
N=fliplr(vander(x))
rref([N y])
function Basis_Int = find_intersection(A, B)
%FIND_INTERSECTION Computes intersection of Col(A) and Col(B)
%   Shows RREF working for manual verification.

fprintf('==============================================\n');
fprintf('   COMPUTING SUBSPACE INTERSECTION\n');
fprintf('==============================================\n');

% --- 1. Setup ---
if size(A, 1) ~= size(B, 1), error('Rows must match.'); end

is_symbolic = isa(A, 'sym') || isa(B, 'sym');
if is_symbolic
    if ~isa(A, 'sym'), A = sym(A); end
    if ~isa(B, 'sym'), B = sym(B); end
    fprintf('Mode: Symbolic (Exact)\n');
else
    fprintf('Mode: Numeric\n');
end

[rows, colsA] = size(A);

% --- 2. Form Matrix M = [A -B] ---
% We want to solve A*c1 = B*c2  =>  A*c1 - B*c2 = 0
if is_symbolic
    M = [A, -B];
else
    M = [A, -B];
end

% ============================================================
%               SHOW WORKING (RREF)
% ============================================================
fprintf('\n--- Working: Solving [A -B]x = 0 ---\n');
fprintf('Matrix [A -B]:\n');
% disp(M); 

if is_symbolic
    R = rref(M);
else
    R = rref(M, 1e-10);
end

fprintf('RREF([A -B]):\n');
if is_symbolic
    disp(R);
else
    disp(rats(R)); % Display as fractions for readability
end
fprintf('Note: Columns without pivots indicate Free Variables.\n');
fprintf('      These free variables generate the intersection.\n');
% ============================================================

% --- 3. Find Null Space & Intersection ---
if is_symbolic
    N = null(M);
else
    N = null(M); 
end

if isempty(N)
    fprintf('\nRESULT: Intersection is {0} (No free variables found).\n');
    Basis_Int = zeros(rows, 1);
    if is_symbolic, Basis_Int = sym(Basis_Int); end
    return;
end

% --- 4. Construct Basis ---
% Take the top half of the null vector (coefficients for A)
% and multiply by A to get the actual vectors in the space.
coeffs_for_A = N(1:colsA, :);
Basis_Int_Raw = A * coeffs_for_A;

% Clean up the basis (remove dependencies)
if is_symbolic
    Basis_Int = simplify(colspace(Basis_Int_Raw));
    % Fallback if colspace doesn't exist in your version: use the raw basis.
    % `norm > 0` would return a sym (undecidable), so test for non-zero
    % entries with isAlways instead.
    if isempty(Basis_Int) && ~all(isAlways(Basis_Int_Raw(:) == 0, 'Unknown', 'false'))
         Basis_Int = Basis_Int_Raw;
    end
else
    Basis_Int = orth(Basis_Int_Raw);
end

fprintf('\n--- Result ---\n');
fprintf('Intersection Dimension: %d\n', size(Basis_Int, 2));
fprintf('Basis for Intersection:\n');
disp(Basis_Int);
fprintf('==============================================\n\n');

end
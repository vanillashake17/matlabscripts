function result = compare_spans(A, B, tol)
%COMPARE_SPANS Compares the spans of column vectors.
%   Handles both Symbolic and Numeric input automatically.
%   Dynamically prints variable names (e.g., "S" and "T") if provided.

% --- 1. Input Validation and Defaults ---
if nargin < 2, error('Requires two input matrices.'); end
if nargin < 3 || isempty(tol), tol = 1e-10; end
if size(A, 1) ~= size(B, 1), error('Rows must match.'); end

% --- 1b. Variable Name Extraction ---
% inputname(n) gets the string name of the nth argument.
nameA = inputname(1);
nameB = inputname(2);

% Fallback: If user passes raw numbers like compare_spans([1 0], [0 1]),
% inputname returns empty. We default to 'A' and 'B' in that case.
if isempty(nameA), nameA = 'A'; end
if isempty(nameB), nameB = 'B'; end

% --- 2. Intelligent Type Handling ---
if isa(A, 'sym') || isa(B, 'sym')
    if ~isa(A, 'sym'), A = sym(A); end
    if ~isa(B, 'sym'), B = sym(B); end
    
    rank_A = rank(A);
    rank_B = rank(B);
    rank_AB = rank([A, B]);
    rank_BA = rank([B, A]);
    
    is_symbolic_mode = true;
else
    rank_A = rank(A, tol);
    rank_B = rank(B, tol);
    rank_AB = rank([A, B], tol);
    rank_BA = rank([B, A], tol);
    
    is_symbolic_mode = false;
end

% --- 3. Determine Relationship ---
is_A_subset_of_B = (rank_BA == rank_B);
is_B_subset_of_A = (rank_AB == rank_A);

if rank_A == 0 && rank_B == 0
    relation_str = 'Trivial zero vector spans';
elseif is_A_subset_of_B && is_B_subset_of_A
    relation_str = 'The spans are equal';
elseif is_A_subset_of_B
    relation_str = sprintf('span(%s) is a proper subset of span(%s)', nameA, nameB);
elseif is_B_subset_of_A
    relation_str = sprintf('span(%s) is a proper subset of span(%s)', nameB, nameA);
else
    relation_str = 'Neither span is a subset of the other';
end

% --- 4. Display Results ---
fprintf('--- Span Relationship Analysis ---\n');
fprintf('Mode: %s\n', string(iff(is_symbolic_mode, "Symbolic (Exact)", "Numeric (Approx)")));
fprintf('Rank(%s): %d, Rank(%s): %d\n', nameA, rank_A, nameB, rank_B);
fprintf('Relation: %s\n\n', relation_str);

% --- Verification Section ---
fprintf('--- Verification ---\n');

if is_symbolic_mode
    rref_AB = rref([A, B]);
    rref_BA = rref([B, A]);
else
    rref_AB = rref([A, B], tol);
    rref_BA = rref([B, A], tol);
end

% Test 1: Is B inside A?
if is_B_subset_of_A
    fprintf('span(%s) IS a subset of span(%s) (Rank [%s %s] = %d matches Rank %s).\n', ...
        nameB, nameA, nameA, nameB, rank_AB, nameA);
else
    fprintf('span(%s) is NOT a subset of span(%s) (Rank [%s %s] = %d != Rank %s).\n', ...
        nameB, nameA, nameA, nameB, rank_AB, nameA);
end
fprintf('RREF([%s %s]):\n', nameA, nameB);
if is_symbolic_mode, disp(rref_AB); else, disp(rats(rref_AB)); end

% Test 2: Is A inside B?
fprintf('\n');
if is_A_subset_of_B
    fprintf('span(%s) IS a subset of span(%s) (Rank [%s %s] = %d matches Rank %s).\n', ...
        nameA, nameB, nameB, nameA, rank_BA, nameB);
else
    fprintf('span(%s) is NOT a subset of span(%s) (Rank [%s %s] = %d != Rank %s).\n', ...
        nameA, nameB, nameB, nameA, rank_BA, nameB);
end
fprintf('RREF([%s %s]):\n', nameB, nameA);
if is_symbolic_mode, disp(rref_BA); else, disp(rats(rref_BA)); end

% Helper
function out = iff(condition, trueVal, falseVal)
    if condition, out = trueVal; else, out = falseVal; end
end

end
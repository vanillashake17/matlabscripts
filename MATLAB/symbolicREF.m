function R = symbolicREF(A)
% symbolicREF Reduces symbolic matrix to row-echelon form WITHOUT destroying singularities.
%   1. Uses Fraction-Free Gaussian Elimination (Cross-Multiplication).
%   2. Simplifies terms using factor(), but NEVER divides the whole row 
%      by a variable (e.g., keeps [a, a] instead of [1, 1]).
%   3. Prioritizes standard pivots (closest to top) to match textbook flow.

    if ~isa(A, 'sym'), A = sym(A); end
    [m, n] = size(A);
    R = A;
    
    pivot_row = 1;
    pivot_col = 1;
    
    while pivot_row <= m && pivot_col <= n
        % --- 1. Conservative Pivot Selection ---
        % Look for the first non-zero candidate. 
        % Do not hunt for "max zeros" as that reshuffles the matrix 
        % too much compared to standard textbook steps.
        
        candidates = pivot_row:m;
        best_row = 0;
        
        % Prefer finding a constant (1, -1) to keep math clean,
        % otherwise just take the first available non-zero row.
        first_nonzero = 0;
        
        for i = candidates
            val = R(i, pivot_col);
            if ~check_is_zero(val)
                if first_nonzero == 0
                    first_nonzero = i;
                end
                
                % If we find a constant, grab it immediately and stop looking
                if isempty(symvar(val))
                    best_row = i;
                    break;
                end
            end
        end
        
        % If no constant found, settle for the first non-zero row
        if best_row == 0
            best_row = first_nonzero;
        end
        
        if best_row == 0
            pivot_col = pivot_col + 1;
            continue;
        end
        
        % --- 2. Swap Rows ---
        if best_row ~= pivot_row
            R([pivot_row, best_row], :) = R([best_row, pivot_row], :);
        end
        
        % --- 3. Fraction-Free Elimination ---
        p = R(pivot_row, pivot_col);
        
        for i = pivot_row + 1 : m
            target = R(i, pivot_col);
            if ~check_is_zero(target)
                % Cross-Multiply: Row = Pivot*Row - Target*PivotRow
                R(i, :) = (p * R(i, :)) - (target * R(pivot_row, :));
                
                % --- CONSERVATIVE SIMPLIFICATION ---
                % Only divide by NUMBERS, never 'a'.
                R(i, :) = conservative_simplify(R(i, :));
            end
        end
        
        pivot_row = pivot_row + 1;
        pivot_col = pivot_col + 1;
    end
end

% -------------------------------------------------------------------------
% THE CONSERVATIVE SIMPLIFIER
% -------------------------------------------------------------------------
function row = conservative_simplify(row)
    % 1. Standard simplify
    row = simplify(row);
    
    % 2. Factor terms (visual cleanup only)
    % This turns a^2+a into a(a+1) without deleting information
    try
        row = factor(row);
    catch
    end
    
    % 3. Numeric GCD Division ONLY
    % We allow dividing [2a, 4a] by 2, but NOT by a.
    idx = find(~check_is_zero(row));
    if isempty(idx), return; end
    vals = row(idx);
    
    try
        % Extract the numeric part of the GCD
        common = vals(1);
        for k = 2:length(vals)
            common = gcd(common, vals(k));
        end
        
        % If the common factor is purely numeric (no variables), divide it out.
        if isempty(symvar(common)) && common ~= 0 && common ~= 1 && common ~= -1
            row = simplify(row / common);
        end
    catch
        % If gcd fails, do nothing
    end
    
    % 4. Sign Normalization (Visual preference)
    first_term = row(idx(1));
    str_val = string(char(expand(first_term)));
    if startsWith(str_val, '-')
        row = -row;
    end
end

function is_z = check_is_zero(val)
    is_z = isequal(simplify(val), sym(0));
end
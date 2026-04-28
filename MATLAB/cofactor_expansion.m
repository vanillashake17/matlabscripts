function [det_val, expansion_str] = cofactor_expansion(A, dim, index)
%COFACTOR_EXPANSION Computes determinant via cofactor expansion.
%   Usage:
%       det_val = cofactor_expansion(A)            % Auto-selects best row/col
%       det_val = cofactor_expansion(A, 'row', 2)  % Expands along Row 2
%       det_val = cofactor_expansion(A, 'col', 1)  % Expands along Col 1
%
%   Returns:
%       det_val: The calculated determinant (symbolic or numeric).
%       expansion_str: A string representing the unsimplified expansion step.

    % --- 1. Input Handling ---
    if ~isa(A, 'sym'), A = sym(A); end
    [n, m] = size(A);
    if n ~= m, error('Matrix must be square.'); end

    % --- 2. Determine Expansion Path ---
    if nargin < 3
        % AUTO MODE: Find row/col with most zeros to minimize work
        row_zeros = sum(A == 0, 2);
        col_zeros = sum(A == 0, 1);
        
        [max_row_z, best_row] = max(row_zeros);
        [max_col_z, best_col] = max(col_zeros);
        
        if max_row_z >= max_col_z
            dim = 'row';
            index = best_row;
            reason = '(Auto-selected: max zeros)';
        else
            dim = 'col';
            index = best_col;
            reason = '(Auto-selected: max zeros)';
        end
    else
        reason = '(User defined)';
    end

    fprintf('==============================================\n');
    fprintf('      COFACTOR EXPANSION (Laplace)            \n');
    fprintf('==============================================\n');
    fprintf('Matrix size: %dx%d\n', n, n);
    fprintf('Expanding along: %s %d %s\n\n', upper(dim), index, reason);

    % --- 3. Compute Expansion ---
    det_val = sym(0);
    terms_str = {};
    
    % Loop through the chosen vector
    for k = 1:n
        if strcmpi(dim, 'row')
            r = index; c = k;
            element = A(r, c);
        else
            r = k; c = index;
            element = A(r, c);
        end
        
        % 3a. Get Sign (-1)^(i+j)
        sign_val = (-1)^(r + c);
        
        % 3b. Get Minor Matrix (Remove row r, col c)
        MinorM = A;
        MinorM(r, :) = [];
        MinorM(:, c) = [];
        
        % 3c. Calculate Determinant of Minor
        det_minor = det(MinorM);
        
        % 3d. Add to total sum
        term_val = sign_val * element * det_minor;
        det_val = det_val + term_val;
        
        % 3e. Formatting for Display
        % Only show terms that aren't strictly zero (unless everything is zero)
        if element ~= 0 || n == 1
            if sign_val == 1
                sign_str = '+';
            else
                sign_str = '-';
            end
            
            % Handle 2x2 minor visualization explicitly
            if size(MinorM, 1) == 2
               min_str = sprintf('|%s %s; %s %s|', string(MinorM(1,1)), string(MinorM(1,2)), ...
                                                   string(MinorM(2,1)), string(MinorM(2,2)));
            elseif size(MinorM, 1) == 1
               min_str = sprintf('(%s)', string(MinorM));
            else
               min_str = sprintf('det(M_%d%d)', r, c);
            end
            
            % Construct string: " + (Element) * |Minor| "
            term_str = sprintf('%s (%s) * %s', sign_str, string(element), min_str);
            terms_str{end+1} = term_str;
            
            % Print step details
            fprintf('Term %d%d: Sign(%s1) * Val(%s) * Det(Minor) = %s\n', ...
                r, c, sign_str, string(element), string(simplify(term_val)));
        end
    end

    % --- 4. Final Output Formatting ---
    full_expansion = strjoin(terms_str, ' ');
    
    % Clean up leading "+" if it exists
    if startsWith(strtrim(full_expansion), '+')
        full_expansion = strtrim(full_expansion);
        full_expansion = full_expansion(2:end);
    end
    
    expansion_str = full_expansion;

    fprintf('\n--- Expansion Expression ---\n');
    disp(expansion_str);
    
    fprintf('\n--- Final Determinant ---\n');
    disp(simplify(det_val));
    fprintf('==============================================\n\n');
end
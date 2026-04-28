function E = orthonormalize_rational(U)
% ORTHONORMALIZE_RATIONAL Converts a set of vectors to an orthonormal set
% using exact symbolic (rational) arithmetic.
%
% Usage: E = orthonormalize_rational(U)
%
% Inputs:
%   U - A matrix where each column is a vector to be orthonormalized.
%
% Outputs:
%   E - A symbolic matrix where each column is an orthonormal vector.
%
% Note: This function will throw an ERROR if the original input vectors
% in U are not strictly orthogonal to each other.

    % Convert input to symbolic to ensure exact rational arithmetic
    U_sym = sym(U);
    [n, k] = size(U_sym);
    E = sym(zeros(n, k));

    % --- Gram-Schmidt Process ---
    for i = 1:k
        % Start with the current vector
        v = U_sym(:, i);
        
        % Subtract projections onto all previously computed orthonormal vectors
        for j = 1:i-1
             v = v - dot(U_sym(:, i), E(:, j)) * E(:, j);
        end
        
        % Check for linear dependence (if v is practically zero)
        if v == sym(zeros(n,1))
             error('Vector u%d is linearly dependent on previous vectors.', i);
        end
        
        % Normalize the resulting vector exactly
        E(:, i) = v / norm(v);
    end
    
    % Display the resulting orthonormal set immediately
    disp('Conversion successful. The orthonormal set is:');
    disp(E);

    % --- Orthogonality Check of ORIGINAL vectors ---
    % We check this after conversion to ensure user sees the result first,
    % but we still strictly enforce the error requirement.
    for i = 1:k
        for j = i+1:k
            d = dot(U_sym(:, i), U_sym(:, j));
            % If dot product is not exactly zero
            if d ~= 0
                fprintf('\nORTHOGONALITY CHECK FAILED:\n');
                fprintf('Original vectors u%d and u%d are NOT orthogonal.\n', i, j);
                fprintf('Dot product = %s\n', char(d));
                error('The original set of vectors is not orthogonal.');
            end
        end
    end
end
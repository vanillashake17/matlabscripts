function result = isOrthogonal(A)
    % Check if the columns of A form an orthogonal set
    % Input: A - matrix where each column is a vector
    
    [~, n] = size(A);
    result = true;
    
    % Check dot product between every pair of columns
    for i = 1:n
        for j = i+1:n
            dot_product = dot(A(:,i), A(:,j));
            if abs(dot_product) > 1e-10  % tolerance for floating point
                fprintf('Columns %d and %d are NOT orthogonal (dot product = %f)\n', i, j, dot_product);
                result = false;
            end
        end
    end
    
    if result
        fprintf('All columns are orthogonal!\n');
    end
end

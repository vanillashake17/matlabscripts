function result = isOrthonormal(A)
    % Check if the columns of A form an orthonormal set
    % (pairwise orthogonal AND each column has unit length),
    % using EXACT symbolic arithmetic so results stay in rational/sqrt form.
    %
    % Input: A - matrix where each column is a vector
    %        (numeric input is auto-converted to sym; pass sym(A) to be explicit)

    if ~isa(A, 'sym')
        A = sym(A);
    end

    [~, n] = size(A);
    result = true;

    % Check unit length of each column (compare ||v||^2 to 1 symbolically)
    fprintf('Column norms:\n');
    for i = 1:n
        nrmSq = simplify(A(:,i).' * A(:,i));
        nrm   = simplify(sqrt(nrmSq));
        fprintf('  ||col %d|| = ', i); disp(nrm);
        if ~isAlways(nrmSq == 1, 'Unknown', 'false')
            fprintf('    -> Column %d is NOT unit length.\n', i);
            result = false;
        end
    end

    % Check dot product between every pair of columns
    fprintf('Pairwise dot products:\n');
    for i = 1:n
        for j = i+1:n
            dp = simplify(A(:,i).' * A(:,j));
            fprintf('  <col %d, col %d> = ', i, j); disp(dp);
            if ~isAlways(dp == 0, 'Unknown', 'false')
                fprintf('    -> Columns %d and %d are NOT orthogonal.\n', i, j);
                result = false;
            end
        end
    end

    if result
        fprintf('All columns are orthonormal (exact).\n');
    else
        fprintf('Columns are NOT orthonormal.\n');
    end
end

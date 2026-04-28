function [x_general, p, Wo, dist] = least_squares(A, v)
% LEAST_SQUARES Performs least-squares analysis showing all RREF workings.
%   Usage: [x, p, Wo, dist] = least_squares(A, b)

    % --- Setup and Validation ---
    % Fix: changed 'symA' to 'sym'
    if ~isa(A, 'sym'), A = sym(A); end
    if ~isa(v, 'sym'), v = sym(v); end
    
    [m, n] = size(A);
    fprintf('==============================================\n');
    fprintf('      LEAST-SQUARES ANALYSIS (With Working)   \n');
    fprintf('==============================================\n\n');

    % --- Step 1: Normal Equations ---
    fprintf('--- Step 1: Formulate Normal Equations ---\n');
    fprintf('Formula: (A''*A)x = (A''*v)\n\n');
    
    A_t = A.';
    A_t_A = A_t * A;
    A_t_v = A_t * v;
    
    fprintf('Matrix (A''*A):\n');
    disp(A_t_A);
    
    fprintf('Vector (A''*v):\n');
    disp(A_t_v);
    
    % --- Step 2: Augmented Matrix & RREF ---
    fprintf('--- Step 2: Solve System [ A''A | A''v ] ---\n');
    
    augmented_matrix = [A_t_A, A_t_v];
    fprintf('Augmented Matrix:\n');
    disp(augmented_matrix);
    
    R = rref(augmented_matrix);
    
    fprintf('RREF of Augmented Matrix:\n');
    disp(R);
    
    % --- Step 3: Analyze Solution ---
    fprintf('--- Step 3: Extract Solution x ---\n');

    if rank(A) == n
        fprintf('Rank matches columns (%d). Unique Solution found.\n', n);
        x_particular = R(:, end);
        x_general = x_particular;
        
        fprintf('\n(1) LEAST-SQUARES SOLUTION [x_hat]:\n');
        disp(x_general);
        fprintf('    (5 s.f.):\n');
        disp(vpa(x_general, 5));
    else
        fprintf('Rank < columns. Infinite Solutions found.\n');
        
        % Find particular solution (xp)
        x_particular = sym(zeros(n, 1));
        pivot_cols = [];
        
        % Identify pivot columns from RREF
        for r = 1:size(R, 1)
            % Find the first non-zero element in the row
            pivot = find(R(r, 1:n), 1, 'first');
            if ~isempty(pivot)
                pivot_cols = [pivot_cols, pivot];
                x_particular(pivot) = R(r, end);
            end
        end
        
        % Find null space basis (homogeneous solutions)
        free_vars = setdiff(1:n, pivot_cols);
        null_basis = [];
        
        for fv_idx = free_vars
            null_vec = sym(zeros(n,1));
            null_vec(fv_idx) = 1; % Set free var to 1
            
            % Back-substitute to find pivot vars
            for r = 1:size(R, 1)
                pivot = find(R(r, 1:n), 1, 'first');
                 if ~isempty(pivot)
                    % pivot_val = -row_val * free_val
                    null_vec(pivot) = -R(r, fv_idx);
                 end
            end
            null_basis = [null_basis, null_vec];
        end
        
        % Construct the general solution string/symbolic object
        syms t [1 length(free_vars)]
        x_general = x_particular;
        
        % Only loop if there are actually free variables
        if ~isempty(free_vars)
            for i = 1:length(free_vars)
                x_general = x_general + t(i)*null_basis(:,i);
            end
        end
        
        fprintf('\n(1) GENERAL LEAST-SQUARES SOLUTION [x_general]:\n');
        printGeneralSolution(x_particular, null_basis, pickParamNames(length(free_vars)));
        fprintf('    (5 s.f.):\n');
        printGeneralSolution(vpa(x_particular, 5), vpa(null_basis, 5), pickParamNames(length(free_vars)));
    end
    
    fprintf('------------------------------------\n\n');
    
    % --- Step 4: Geometric Results ---
    % Use only the particular solution for geometric projection
    p = simplify(A * x_particular);
    Wo = simplify(v - p);
    dist = simplify(norm(Wo));

    fprintf('--- Step 4: Geometric Projections ---\n');
    fprintf('Projection p = A * x_particular\n');
    fprintf('\n(2) ORTHOGONAL PROJECTION [p]:\n');
    disp(p);
    fprintf('    (5 s.f.):\n');
    disp(vpa(p, 5));

    fprintf('(3) ORTHOGONAL COMPONENT [v - p]:\n');
    disp(Wo);
    fprintf('    (5 s.f.):\n');
    disp(vpa(Wo, 5));

    fprintf('(4) LEAST-SQUARES ERROR (Distance):\n');
    disp(dist);
    fprintf('    (5 s.f.):\n');
    disp(vpa(dist, 5));
    fprintf('==============================================\n');
end


function names = pickParamNames(k)
    if k == 0
        names = {};
    elseif k <= 3
        all_names = {'s', 't', 'r'};
        names = all_names(1:k);
    else
        names = arrayfun(@(i) sprintf('t%d', i), 1:k, 'UniformOutput', false);
    end
end


function printGeneralSolution(xp, basis, paramNames)
    xp = xp(:);
    n = length(xp);
    k = size(basis, 2);

    xpStr = cell(n, 1);
    for i = 1:n
        xpStr{i} = strtrim(char(xp(i)));
    end

    basisStr = cell(n, k);
    for j = 1:k
        for i = 1:n
            basisStr{i, j} = strtrim(char(basis(i, j)));
        end
    end

    wxp = max(cellfun(@length, xpStr));
    wb = zeros(1, k);
    for j = 1:k
        wb(j) = max(cellfun(@length, basisStr(:, j)));
    end

    if k > 0
        wp = max(cellfun(@length, paramNames));
    else
        wp = 0;
    end

    labelRow = max(1, ceil(n / 2));

    fprintf('\n');
    for i = 1:n
        fprintf('  [ %*s ]', wxp, xpStr{i});
        for j = 1:k
            if i == labelRow
                fprintf(' + %-*s', wp, paramNames{j});
            else
                fprintf('   %-*s', wp, '');
            end
            fprintf(' [ %*s ]', wb(j), basisStr{i, j});
        end
        fprintf('\n');
    end

    if k > 0
        fprintf('  where %s in R\n', strjoin(paramNames, ', '));
    end
    fprintf('\n');
end
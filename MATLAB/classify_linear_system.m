function cases = classify_linear_system(varargin)
% CLASSIFY_LINEAR_SYSTEM  Classifies a parametric linear system using a
% robust, case-by-case methodology.
%
%   classify_linear_system(M)      where M = [A | b] is the augmented matrix
%                                  (any m-by-(n+1) shape — square, over-,
%                                   or under-determined all OK)
%   classify_linear_system(A, b)   pass A and b separately, like Ax = b
%
% To classify the HOMOGENEOUS system Ax = 0 for a parametric coefficient
% matrix A, call:  classify_linear_system(A, zeros(size(A,1), 1))

    if nargin == 1
        M = sym(varargin{1});
        % Sanity guard against the easy mistake of passing a coefficient
        % matrix and expecting Ax = 0 to be analysed automatically.
        [r, c] = size(M);
        if r == c
            warning('classify_linear_system:ambiguousInput', ...
                ['Input is %dx%d (square). It is being treated as an ', ...
                 'AUGMENTED matrix [A | b] with %d unknowns and %d equations. ', ...
                 'For the homogeneous system Ax=0, call ', ...
                 'classify_linear_system(A, zeros(%d, 1)) instead.'], ...
                r, c, c-1, r, r);
        end
    elseif nargin == 2
        A = sym(varargin{1});
        b = sym(varargin{2});
        b = b(:);
        assert(size(A, 1) == size(b, 1), ...
            'A and b must have the same number of rows.');
        M = [A, b];
    else
        error('Usage: classify_linear_system(M) or classify_linear_system(A, b)');
    end

    % --- 1. SETUP AND VALIDATION ---
    [num_eqs, num_cols] = size(M);
    n = num_cols - 1; % Number of variables (coefficient columns of M)
    if num_cols < 2
        error(['Input must have at least 2 columns: at least one ', ...
               'coefficient column plus the RHS. Got %d-by-%d.'], ...
              num_eqs, num_cols);
    end

    % Find symbolic parameters in the matrix.
    vars = symvar(M);
    if isempty(vars)
        fprintf('--- System has no symbolic variables. Performing standard analysis. ---\n');
        cases.no_symbolic_vars = analyze_case(M, 'Constant System');
        return;
    end
    if numel(vars) > 4
        error(['This function supports up to 4 symbolic parameters. ', ...
               'Got %d (%s).'], numel(vars), char(join(string(vars), ', ')));
    end

    multi_param = numel(vars) > 1;
    p = vars(1); % Primary parameter (used for headings + 1-param fast path).

    if multi_param
        fprintf('--- Classifying Linear System for Parameters %s ---\n\n', ...
                char(join(string(vars), ', ')));
    else
        fprintf('--- Classifying Linear System for Parameter ''%s'' ---\n\n', char(p));
    end
    fprintf('Input Augmented Matrix M:\n');
    disp(M);

    % --- 2. FIND CRITICAL CONDITIONS ---
    % A "critical substitution" is a (var, value) pair such that setting
    % var := value zeros out a critical condition (denominator, consistency
    % row, or det(A)). Each substitution becomes one special case. For 1
    % parameter this collapses to the classical "list of critical values".
    fprintf('Step 1: Finding critical conditions by inspecting the general RREF...\n');
    R_general_raw = rref(M);

    crit_subs = struct('var', {}, 'val', {});   % each special case to enumerate
    crit_conditions = sym([]);                   % printed as "expr = 0"

    % (a) Denominators of the parametric RREF — values where pivot
    %     assumptions break (rank can change).
    [~, d] = numden(R_general_raw);
    denominators = unique(d(d~=1));
    if ~isempty(denominators)
        fprintf('Found denominators in the general RREF: %s\n', ...
                join(string(denominators), ', '));
        for i = 1:numel(denominators)
            crit_conditions = [crit_conditions; denominators(i)]; %#ok<AGROW>
            crit_subs = [crit_subs, expand_critical(denominators(i), vars)]; %#ok<AGROW>
        end
    end

    % (b) Consistency rows: [0 ... 0 | f(params)]. Consistent iff f = 0,
    %     so each root of f is a boundary between consistent / inconsistent.
    n_eqs = size(R_general_raw, 1);
    A_part = R_general_raw(:, 1:n);
    b_part = R_general_raw(:, n+1);
    consistency_exprs = sym([]);
    for i = 1:n_eqs
        row_zero = all(arrayfun( ...
            @(x) isAlways(x == 0, 'Unknown', 'false'), A_part(i, :)));
        bi = simplify(b_part(i));
        if row_zero && ~isAlways(bi == 0, 'Unknown', 'false') ...
                    && ~isempty(symvar(bi))
            consistency_exprs = [consistency_exprs; bi]; %#ok<AGROW>
            crit_conditions = [crit_conditions; bi]; %#ok<AGROW>
            crit_subs = [crit_subs, expand_critical(bi, vars)]; %#ok<AGROW>
        end
    end
    if ~isempty(consistency_exprs)
        fprintf('Consistency-row expressions in the RREF: %s\n', ...
                join(string(unique(consistency_exprs)), ', '));
    end

    % (c) Rank-drop conditions for the coefficient matrix.
    %     Square A:        det(A) = 0 marks the rank-drop boundary.
    %     Tall A  (m > n): det(A'*A) = sum of squared n-by-n minors
    %                      (Cauchy–Binet); zero iff rank(A) < n.
    %     Wide A  (m < n): det(A*A') is the same idea for full row rank;
    %                      unique solutions are impossible regardless, but
    %                      this still marks where the consistency regime
    %                      changes.
    %     Without this branch for non-square A, RREF can silently divide
    %     by parametric pivots whose vanishing isn't visible in the final
    %     denominators (e.g. AY2425 makeup midterm Q1: 5x4 matrix where
    %     a=0 and a=b drop the rank but classify only saw a-1 from the
    %     RREF denominators).
    A_only = M(:, 1:n);
    [m_rows, n_cols] = size(A_only);
    if m_rows == n_cols
        rank_poly = simplify(det(A_only));
        rank_label = 'det(A)';
    elseif m_rows > n_cols
        rank_poly = simplify(det(A_only.' * A_only));
        rank_label = 'det(A''*A)';
    else
        rank_poly = simplify(det(A_only * A_only.'));
        rank_label = 'det(A*A'')';
    end
    if ~isAlways(rank_poly == 0, 'Unknown', 'false') && ~isempty(symvar(rank_poly))
        fprintf('%s = %s\n', rank_label, char(rank_poly));
        crit_conditions = [crit_conditions; rank_poly];
        crit_subs = [crit_subs, expand_critical(rank_poly, vars)];
    elseif isAlways(rank_poly == 0, 'Unknown', 'false')
        % Fallback: the primary rank polynomial is IDENTICALLY zero, meaning
        % A is generically rank-deficient (e.g. a column is zero for every
        % value of the parameters). Rank can still drop FURTHER on a
        % parameter variety, which would push the system from "infinite
        % solutions with k parameters" to "infinite solutions with k+1
        % parameters". Enumerate all r×r minors (where r is the generic
        % rank) and treat each non-trivial one as a rank-drop boundary.
        % Without this, e.g. the 4x4 system from AY2425 final Q1 misses the
        % a=0,b=0 joint case that yields 2-parameter infinite solutions.
        try
            r_gen = double(rank(A_only));
        catch
            r_gen = min(m_rows, n_cols) - 1;
        end
        if r_gen > 0 && r_gen < min(m_rows, n_cols)
            row_combos = nchoosek(1:m_rows, r_gen);
            col_combos = nchoosek(1:n_cols, r_gen);
            minor_polys = sym([]);
            for ri = 1:size(row_combos, 1)
                for ci = 1:size(col_combos, 1)
                    sub = A_only(row_combos(ri, :), col_combos(ci, :));
                    d = simplify(det(sub));
                    if ~isAlways(d == 0, 'Unknown', 'false') ...
                            && ~isempty(symvar(d))
                        minor_polys = [minor_polys; d]; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(minor_polys)
                minor_polys = unique(minor_polys);
                fprintf('%dx%d minors with parametric rank-drop: %s\n', ...
                    r_gen, r_gen, ...
                    char(join(string(minor_polys), ', ')));
                for k = 1:numel(minor_polys)
                    crit_conditions = [crit_conditions; minor_polys(k)]; %#ok<AGROW>
                    crit_subs = [crit_subs, ...
                        expand_critical(minor_polys(k), vars)]; %#ok<AGROW>
                end
            end
        end
    end

    % (d) Left null space of the ORIGINAL coefficient matrix. Any vector
    %     w with w'*A = 0 forces a consistency requirement w'*b = 0.
    %     This catches conditions that MATLAB's symbolic rref hides when
    %     it normalizes a parametric divisor (e.g. turning [0...0 | b-2]
    %     into [0...0 | 1] and erasing the b=2 condition).
    b_part_orig = M(:, n+1);
    try
        W = null(A_only.');
    catch
        W = sym([]);
    end
    if ~isempty(W)
        ln_conds = sym([]);
        for k = 1:size(W, 2)
            cond = simplify(W(:, k).' * b_part_orig);
            if ~isAlways(cond == 0, 'Unknown', 'false') ...
                       && ~isempty(symvar(cond))
                ln_conds = [ln_conds; cond]; %#ok<AGROW>
                crit_conditions = [crit_conditions; cond]; %#ok<AGROW>
                crit_subs = [crit_subs, expand_critical(cond, vars)]; %#ok<AGROW>
            end
        end
        if ~isempty(ln_conds)
            fprintf('Left-null consistency conditions: %s\n', ...
                join(string(unique(ln_conds)), ', '));
        end
    end

    % Deduplicate substitutions (same (var, value) pair from multiple sources).
    crit_subs = dedupe_subs(crit_subs);

    % For backward compatibility with the 1-param case, also expose
    % critical_values as a vector of values of the primary parameter.
    critical_values = sym([]);
    for k = 1:numel(crit_subs)
        if isequal(crit_subs(k).var, p)
            critical_values = [critical_values; crit_subs(k).val]; %#ok<AGROW>
        end
    end
    
    if isempty(crit_subs)
        fprintf('No critical conditions found. Only a general case analysis is needed.\n\n');
        general_cond_str = 'General Case (no special cases exist)';
    else
        crit_conditions_uniq = unique(simplify(crit_conditions));
        sub_strs  = arrayfun(@(s) string( ...
            sprintf('%s = %s', char(crit_subs(s).var), short_string(crit_subs(s).val))), ...
            1:numel(crit_subs));
        fprintf('Critical conditions (system behavior may change when any of these vanish):\n');
        for c = 1:numel(crit_conditions_uniq)
            fprintf('   %s = 0\n', char(crit_conditions_uniq(c)));
        end
        fprintf('Enumerated as substitutions: %s\n\n', join(sub_strs, ', '));
        general_cond_str = 'General Case (none of the critical conditions vanish)';
    end

    % --- 3. ANALYZE THE GENERAL CASE ---
    cases.general_case = analyze_case(R_general_raw, general_cond_str);

    % Track all leaf cases for the rolled-up summary table.
    % Each leaf has: .conditions (cell of strings) and .classification.
    cases.leaves = {};
    if isempty(crit_subs)
        cases.leaves{end+1} = struct( ...
            'conditions', {{'(no parameters)'}}, ...
            'classification', cases.general_case.classification);
    else
        gen_neg_strs = arrayfun(@(s) ...
            sprintf('%s ≠ %s', char(crit_subs(s).var), short_string(crit_subs(s).val)), ...
            1:numel(crit_subs), 'UniformOutput', false);
        cases.leaves{end+1} = struct( ...
            'conditions', {gen_neg_strs}, ...
            'classification', cases.general_case.classification);
    end

    % --- 4. ANALYZE EACH SPECIAL CASE -------------------------------------
    % Two flavours:
    %   (a) atomic — substitution leaves no remaining parameters; classify
    %       directly via analyze_case and bucket by NO/INFINITE/UNIQUE.
    %   (b) parametric remainder — substitution still has params; recurse via
    %       classify_linear_system to enumerate the JOINT cases (e.g. a=0
    %       AND b=2 in addition to a=0 alone). Recursive output is shown
    %       verbatim, not bucketed.
    cases.special_cases = struct();
    if ~isempty(crit_subs)
        buckets = struct( ...
            'NoSolution',         {{}}, ...
            'InfiniteSolutions1', {{}}, ...
            'InfiniteSolutions2', {{}}, ...
            'InfiniteSolutionsK', {{}}, ...
            'UniqueSolution',     {{}});
        recursive_entries = {};   % cell array of {label, fieldname, output, result}

        for i = 1:numel(crit_subs)
            sv = crit_subs(i);
            v_str = short_string(sv.val);
            sub_label = sprintf('%s = %s', char(sv.var), v_str);
            safe_fieldname = matlab.lang.makeValidName( ...
                sprintf('case_for_%s_eq_%s', char(sv.var), v_str));

            M_special = subs(M, sv.var, sv.val);

            if isempty(symvar(M_special))
                % Atomic case — substitution removed the last parameter.
                special_cond_str = sprintf('Special Case (%s)', sub_label); %#ok<NASGU>
                M_special_local = M_special; %#ok<NASGU>
                output_str = evalc( ...
                    'r = analyze_case(M_special_local, special_cond_str);');
                if strcmp(r.classification, 'No Solution')
                    key = 'NoSolution';
                elseif strcmp(r.classification, 'Unique Solution')
                    key = 'UniqueSolution';
                elseif startsWith(r.classification, 'Infinite Solutions')
                    if r.num_params == 1
                        key = 'InfiniteSolutions1';
                    elseif r.num_params == 2
                        key = 'InfiniteSolutions2';
                    else
                        key = 'InfiniteSolutionsK';
                    end
                else
                    key = 'UniqueSolution';
                end
                buckets.(key){end+1} = struct( ...
                    'output', output_str, ...
                    'fieldname', safe_fieldname, ...
                    'sub', sv, ...
                    'label', sub_label, ...
                    'result', r);
            else
                % Parametric remainder — recurse on the substituted matrix.
                % This catches joint conditions like (a=0 AND b=2).
                output_str = evalc( ...
                    ['fprintf(''\n========================================' ...
                     '\n  Special Case (%s)  — still parametric, recursing' ...
                     '\n========================================\n'', sub_label);' ...
                     'r = classify_linear_system(M_special);']);
                recursive_entries{end+1} = struct( ...
                    'output', output_str, ...
                    'fieldname', safe_fieldname, ...
                    'sub', sv, ...
                    'label', sub_label, ...
                    'result', r); %#ok<AGROW>
            end
        end

        group_order = {'NoSolution', 'InfiniteSolutions1', ...
                       'InfiniteSolutions2', 'InfiniteSolutionsK', ...
                       'UniqueSolution'};
        group_label = struct( ...
            'NoSolution',         'NO SOLUTION', ...
            'InfiniteSolutions1', 'INFINITE SOLUTIONS (1 PARAMETER s)', ...
            'InfiniteSolutions2', 'INFINITE SOLUTIONS (2 PARAMETERS s, t)', ...
            'InfiniteSolutionsK', 'INFINITE SOLUTIONS (3+ PARAMETERS)', ...
            'UniqueSolution',     'UNIQUE SOLUTION');

        fprintf('\n--- Special Case Summary ---\n');
        for g = 1:numel(group_order)
            key = group_order{g};
            entries = buckets.(key);
            if isempty(entries), continue; end
            labels = arrayfun(@(k) string(entries{k}.label), 1:numel(entries));
            fprintf('  %-20s %s\n', [group_label.(key) ':'], join(labels, ',   '));
        end
        if ~isempty(recursive_entries)
            labels = arrayfun(@(k) string(recursive_entries{k}.label), ...
                              1:numel(recursive_entries));
            fprintf('  %-20s %s\n', 'STILL PARAMETRIC:', join(labels, ',   '));
            fprintf('   (each is analyzed recursively below — joint conditions handled there)\n');
        end

        for g = 1:numel(group_order)
            key = group_order{g};
            entries = buckets.(key);
            if isempty(entries), continue; end
            fprintf('\n========== %s cases ==========\n', group_label.(key));
            for k = 1:numel(entries)
                fprintf('%s', entries{k}.output);
                cases.special_cases.(entries{k}.fieldname) = entries{k}.result;
                % Atomic leaf: just the substitution as condition.
                sv = entries{k}.sub;
                cases.leaves{end+1} = struct( ...
                    'conditions', {{sprintf('%s = %s', char(sv.var), short_string(sv.val))}}, ...
                    'classification', entries{k}.result.classification);
            end
        end
        for k = 1:numel(recursive_entries)
            fprintf('%s', recursive_entries{k}.output);
            cases.special_cases.(recursive_entries{k}.fieldname) = recursive_entries{k}.result;
            % Recursive leaves: prepend my substitution to each inner leaf.
            sv = recursive_entries{k}.sub;
            my_cond = sprintf('%s = %s', char(sv.var), short_string(sv.val));
            inner = recursive_entries{k}.result.leaves;
            for j = 1:numel(inner)
                merged = inner{j};
                % If inner condition is "(no parameters)", replace with my_cond.
                if numel(merged.conditions) == 1 && ...
                        strcmp(merged.conditions{1}, '(no parameters)')
                    merged.conditions = {my_cond};
                else
                    merged.conditions = [{my_cond}, merged.conditions];
                end
                cases.leaves{end+1} = merged;
            end
        end
    end

    % --- 5. ROLL-UP SUMMARY TABLE (top-level call only) -------------------
    stack = dbstack;
    is_recursive = numel(stack) >= 2 && ...
                   strcmp(stack(2).name, 'classify_linear_system');
    if ~is_recursive
        cases.leaves = dedupe_leaves(cases.leaves);
        print_leaves_table(cases.leaves);
    end

    fprintf('\n--- Classification Complete ---\n');
end

% -------------------------------------------------------------------------
function subs_list = expand_critical(eq, vars)
% Convert a critical-condition expression `eq` (a polynomial in some subset
% of `vars`) into a list of (var, value) substitution structs that zero it
% out. Strategy: factor `eq` first; for each factor that depends on at least
% one variable in `vars`, solve for one of those variables.
%
% For a single-variable factor like (a - 1)(a + 2) -> [(a, 1), (a, -2)].
% For a multi-variable factor like (a*b - 1) -> picks one variable and
% solves: solve(a*b - 1, a) = 1/b, giving substitution (a, 1/b).
    subs_list = struct('var', {}, 'val', {});
    if isAlways(eq == 0, 'Unknown', 'false') || isempty(symvar(eq))
        return;
    end

    % factor() returns either a single sym (already irreducible) or a vector.
    try
        fs = factor(simplify(eq));
    catch
        fs = simplify(eq);
    end
    if numel(fs) == 0
        fs = simplify(eq);
    end

    for i = 1:numel(fs)
        f = fs(i);
        fvars = symvar(f);
        if isempty(fvars)
            continue; % constant factor -> nothing to enumerate
        end
        % Pick the variable to solve for: prefer ones with lowest degree in f
        % (often gives a clean closed-form root).
        chosen = fvars(1);
        try
            best_deg = polynomialDegree(f, fvars(1));
        catch
            best_deg = Inf;
        end
        for j = 2:numel(fvars)
            try
                d_here = polynomialDegree(f, fvars(j));
                if d_here < best_deg
                    chosen = fvars(j);
                    best_deg = d_here;
                end
            catch
            end
        end

        try
            sols = solve(f == 0, chosen);
        catch
            continue;
        end
        if isempty(sols), continue; end

        % Filter to real solutions only. For purely numeric solutions, skip
        % imaginary roots. For symbolic-in-other-vars solutions, skip those
        % that carry the imaginary unit (e.g. det(A'*A) factors like
        % 12a^2 + 4ab + 11b^2 have negative discriminant and produce
        % a = -b/6 ± (2√2/3)i·b — never real for real (a,b) except (0,0),
        % which is already reached via other factors).
        for k = 1:numel(sols)
            v = sols(k);
            if isempty(symvar(v))
                try
                    if imag(vpa(v, 12)) ~= 0
                        continue;
                    end
                catch
                end
            else
                try
                    if has(v, sym(1i))
                        continue;
                    end
                catch
                end
            end
            subs_list(end+1) = struct('var', chosen, 'val', v); %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function out = dedupe_subs(subs_list)
% Collapse duplicate (var, value) substitutions. Two subs are duplicates if
% they share the same variable and their values are equal symbolically.
    out = struct('var', {}, 'val', {});
    for i = 1:numel(subs_list)
        is_dup = false;
        for j = 1:numel(out)
            if isequal(subs_list(i).var, out(j).var) && ...
               isAlways(subs_list(i).val == out(j).val, 'Unknown', 'false')
                is_dup = true;
                break;
            end
        end
        if ~is_dup
            out(end+1) = subs_list(i); %#ok<AGROW>
        end
    end
end

function result = analyze_case(Matrix, case_name)
% Helper function to analyze a specific matrix and print a report.
    fprintf('\n--- %s ---\n', case_name);
    if ~isequal(Matrix, rref(Matrix))
        fprintf('Augmented Matrix for this case:\n'); disp(Matrix);
        R = rref(Matrix);
        fprintf('\nRREF for this case:\n'); disp(R);
    else
        R = Matrix; fprintf('RREF for this case:\n'); disp(R);
    end
    n = size(R, 2) - 1; A_part = R(:, 1:n); b_part = R(:, n+1);
    is_consistent = true;
    for i = 1:size(R, 1)
        if all(A_part(i,:) == 0) && b_part(i) ~= 0, is_consistent = false; break; end
    end
    if ~is_consistent
        result.classification = 'No Solution';
        fprintf('Analysis: A row of the form [0 0 ... | non-zero] was found.\n');
        fprintf('Conclusion: The system is inconsistent and has NO SOLUTION.\n');
    else
        % --- FIX: Manually find pivot columns for compatibility ---
        pivot_cols = [];
        for i = 1:size(A_part, 1)
            % Find the column index of the first non-zero element in each row
            col_idx = find(A_part(i,:), 1); 
            if ~isempty(col_idx)
                pivot_cols(end+1) = col_idx;
            end
        end
        rankA = length(pivot_cols);
        % --- END OF FIX ---
        
        if rankA == n
            result.classification = 'Unique Solution';
            fprintf('Analysis: The system is consistent and there are no free variables (rank = %d = n).\n', n);
            fprintf('Conclusion: The system has a UNIQUE SOLUTION.\n');
            result.solution = simplify(R(1:n, n+1));
            fprintf('Solution:\n'); disp(result.solution);
        else
            num_params = n - rankA;
            result.num_params = num_params;
            free_cols = setdiff(1:n, pivot_cols);
            param_names = pick_param_names(numel(free_cols));
            if num_params == 1
                result.classification = 'Infinite Solutions (1 parameter)';
            else
                result.classification = sprintf( ...
                    'Infinite Solutions (%d parameters)', num_params);
            end
            fprintf('Analysis: The system is consistent and has %d free variable(s) (rank = %d < n).\n', num_params, rankA);
            fprintf('Conclusion: The system has INFINITE SOLUTIONS (%d parameter(s): %s).\n', ...
                num_params, strjoin(param_names, ', '));
            params = sym(param_names);

            % Build the general solution: free vars take parameter symbols;
            % pivot vars are expressed via their RREF row.
            x_gen = sym(zeros(n, 1));
            for k = 1:numel(free_cols)
                x_gen(free_cols(k)) = params(k);
            end
            for r = 1:rankA
                pc = pivot_cols(r);
                val = R(r, n+1);
                for k = 1:numel(free_cols)
                    val = val - R(r, free_cols(k)) * params(k);
                end
                x_gen(pc) = simplify(val);
            end

            % Decompose into particular + null-space basis.
            x_part = subs(x_gen, params, sym(zeros(1, numel(params))));
            null_basis = sym(zeros(n, numel(params)));
            for k = 1:numel(params)
                e_k = sym(zeros(1, numel(params)));
                e_k(k) = 1;
                null_basis(:, k) = subs(x_gen, params, e_k) - x_part;
            end

            free_var_labels = arrayfun( ...
                @(k) sprintf('x_%d = %s', free_cols(k), param_names{k}), ...
                1:numel(free_cols), 'UniformOutput', false);
            fprintf('Free parameter(s): %s   (%s)\n', ...
                strjoin(param_names, ', '), strjoin(free_var_labels, ', '));
            fprintf('General solution x = x_p + %s:\n', ...
                strjoin(arrayfun(@(k) sprintf('%s*v_%d', param_names{k}, k), ...
                                  1:numel(params), 'UniformOutput', false), ' + '));
            disp(simplify(x_gen));
            fprintf('Particular solution x_p (all parameters = 0):\n');
            disp(x_part);
            fprintf('Null-space basis (columns):\n');
            disp(null_basis);

            result.solution     = simplify(x_gen);
            result.particular   = x_part;
            result.null_basis   = null_basis;
            result.parameters   = params;
        end
    end
    result.RREF = R;
end


function names = pick_param_names(k)
    if k == 0
        names = {};
    elseif k <= 3
        all_names = {'s', 't', 'r'};
        names = all_names(1:k);
    else
        names = arrayfun(@(i) sprintf('t%d', i), 1:k, 'UniformOutput', false);
    end
end

function s = short_string(v)
% Return a short, human-friendly string for a symbolic scalar.
% Strategy: prefer the symbolic form (e.g. "sqrt(2)", "1/3", "-2") when it is
% reasonably short. Otherwise fall back to a 6-decimal-place numeric form.
    s = char(v);
    if numel(s) <= 20
        return;
    end
    try
        d = double(v);
        if isfinite(d)
            s = sprintf('%.6f', d);
            % Trim trailing zeros for clean display: 1.414214 stays;
            % 2.000000 becomes 2; 1.500000 becomes 1.5.
            if contains(s, '.')
                s = regexprep(s, '0+$', '');
                s = regexprep(s, '\.$', '');
            end
        end
    catch
        % If we cannot get a numeric value, leave the long symbolic form.
    end
end

% -------------------------------------------------------------------------
function out = dedupe_leaves(leaves)
% Collapse leaves with identical condition sets (order-independent) and
% identical classification. Different recursion paths can both reach the
% same joint condition (e.g. (a=0, b=2) is reached via a=0→b=2 AND b=2→a=0).
    out = {};
    for i = 1:numel(leaves)
        L = leaves{i};
        key_i = canonical_key(L);
        is_dup = false;
        for j = 1:numel(out)
            if strcmp(canonical_key(out{j}), key_i)
                is_dup = true;
                break;
            end
        end
        if ~is_dup
            out{end+1} = L; %#ok<AGROW>
        end
    end
end

function k = canonical_key(L)
% Sort conditions to make key order-independent, then concat with the class.
    conds_sorted = sort(string(L.conditions));
    k = char(strjoin(conds_sorted, ' & ') + " || " + string(L.classification));
end

% -------------------------------------------------------------------------
function print_leaves_table(leaves)
% Print a clean rolled-up table of all (joint condition -> classification)
% leaves found anywhere in the recursion tree.
    if isempty(leaves), return; end

    % Sort each leaf's conditions alphabetically so equivalent leaves print
    % in the same order regardless of which recursion path produced them.
    cond_strs = cell(1, numel(leaves));
    cond_w = strlength('Conditions');
    for i = 1:numel(leaves)
        sorted = sort(string(leaves{i}.conditions));
        cond_strs{i} = char(strjoin(sorted, ',  '));
        if strlength(cond_strs{i}) > cond_w
            cond_w = strlength(cond_strs{i});
        end
    end

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('  CLASSIFICATION SUMMARY (joint conditions, leaf cases)\n');
    fprintf('============================================================\n');
    fprintf('  %-5s  %-*s  %s\n', 'Case', cond_w, 'Conditions', 'Classification');
    fprintf('  %s  %s  %s\n', repmat('-', 1, 5), ...
        repmat('-', 1, cond_w), repmat('-', 1, 20));

    % Build display order dynamically: No Solution → Infinite Solutions
    % (sorted by parameter count, ascending) → Unique Solution. Each
    % Infinite Solutions classification carries its own parameter count
    % from analyze_case (e.g. 'Infinite Solutions (1 parameter)'); collect
    % the unique strings present in the leaves and sort them numerically.
    inf_cls = {};
    for i = 1:numel(leaves)
        c = leaves{i}.classification;
        if startsWith(c, 'Infinite Solutions')
            inf_cls{end+1} = c; %#ok<AGROW>
        end
    end
    inf_cls = unique(inf_cls);
    if ~isempty(inf_cls)
        nums = zeros(1, numel(inf_cls));
        for i = 1:numel(inf_cls)
            tok = regexp(inf_cls{i}, '\((\d+)\s*parameter', 'tokens', 'once');
            if isempty(tok)
                nums(i) = inf;
            else
                nums(i) = str2double(tok{1});
            end
        end
        [~, idx] = sort(nums);
        inf_cls = inf_cls(idx);
    end
    order = [{'No Solution'}, inf_cls, {'Unique Solution'}];
    label_map = containers.Map();
    label_map('No Solution') = 'NO SOLUTION';
    label_map('Unique Solution') = 'UNIQUE SOLUTION';
    for i = 1:numel(inf_cls)
        label_map(inf_cls{i}) = upper(inf_cls{i});
    end

    n = 1;
    for g = 1:numel(order)
        cls = order{g};
        for i = 1:numel(leaves)
            if strcmp(leaves{i}.classification, cls)
                fprintf('  (%-3d)  %-*s  %s\n', n, cond_w, cond_strs{i}, label_map(cls));
                n = n + 1;
            end
        end
    end
    fprintf('============================================================\n');
end
function x = solveLinearODE(A, x0, t0)
    % Solve y' = A y for constant matrix A using the Jordan form.
    % Walks Jordan blocks manually so output is real for complex eigenvalues
    % (no leftover 1i) and chain solutions are emitted in textbook form.
    %
    %   solveLinearODE(A)             general solution with arbitrary constants c1..cn
    %   solveLinearODE(A, x0)         IVP with x(0) = x0
    %   solveLinearODE(A, x0, t0)     IVP with x(t0) = x0
    %
    % Prints each fundamental solution x_k(t) individually, then the
    % general or particular solution. Returns the solution if requested.

    if nargin < 2, x0 = []; end
    if nargin < 3, t0 = 0; end

    syms t real
    A = sym(A);
    n = size(A, 1);

    [V, J] = jordan(A);

    % --- Identify Jordan block ranges ---------------------------------
    block_starts = 1;
    for i = 2:n
        if J(i-1, i) == 0
            block_starts(end+1) = i; %#ok<AGROW>
        end
    end
    block_ends = [block_starts(2:end) - 1, n];
    nblocks = numel(block_starts);

    block_lambdas = sym(zeros(1, nblocks));
    for b = 1:nblocks
        block_lambdas(b) = J(block_starts(b), block_starts(b));
    end

    % --- Build fundamental solutions block by block -------------------
    used = false(1, nblocks);
    fundamentals = sym(zeros(n, 0));

    for b = 1:nblocks
        if used(b), continue; end
        used(b) = true;

        bs = block_starts(b);
        bsize = block_ends(b) - bs + 1;
        lam = block_lambdas(b);

        im_lam = simplify(imag(lam));
        is_real = isAlways(im_lam == 0, 'Unknown', 'false');

        if is_real
            % Chain solutions: x_j(t) = e^{lam t} * sum_{i=0}^{j-1} (t^i/i!) v_{j-i}
            chain = build_real_chain(V, bs, bsize, lam, t);
            fundamentals = [fundamentals, chain]; %#ok<AGROW>
        else
            % Pair with conjugate block, emit two real chains
            cj_b = find_conjugate_block(block_lambdas, used, lam);
            if isempty(cj_b)
                % Defensive fallback: emit raw complex chain
                fundamentals = [fundamentals, build_real_chain(V, bs, bsize, lam, t)]; %#ok<AGROW>
                continue;
            end
            used(cj_b) = true;

            alpha = simplify(real(lam));
            beta  = simplify(imag(lam));

            [reals, imags] = build_complex_pair(V, bs, bsize, alpha, beta, t);
            fundamentals = [fundamentals, reals, imags]; %#ok<AGROW>
        end
    end

    % --- Clear rational denominators per fundamental column -----------
    % jordan(A) sometimes returns eigenvectors with fractional entries
    % (e.g. [1/2; 1/2; 1]). Each fundamental solution is unique up to a
    % nonzero scalar, so multiply through by the LCM of constant rational
    % denominators to keep workings in integer form.
    for k = 1:size(fundamentals, 2)
        fundamentals(:, k) = clear_rational_fractions(fundamentals(:, k), t);
    end

    % --- Print fundamental solutions ----------------------------------
    fprintf('Fundamental solutions:\n');
    for k = 1:size(fundamentals, 2)
        fprintf('\n  x_%d(t) =\n', k);
        disp(fundamentals(:, k));
    end

    % --- General or particular solution -------------------------------
    nfund = size(fundamentals, 2);
    if isempty(x0)
        c = sym('c', [nfund 1]);
        % No simplify on the matrix-vector product — it would fold trig sums
        % into phase form. Symbolic auto-eval gives the textbook expansion.
        xlocal = fundamentals * c;
        fprintf('\nGeneral solution x(t) = c_1*x_1(t) + ... + c_%d*x_%d(t):\n', ...
                nfund, nfund);
        disp(xlocal);
    else
        x0v = sym(x0(:));
        Phi_t0 = simplify(subs(fundamentals, t, sym(t0)));
        c_vals = simplify(Phi_t0 \ x0v);
        fprintf('\nConstants from x(%s) = %s:\n', char(sym(t0)), char(x0v.'));
        for k = 1:numel(c_vals)
            fprintf('  c_%d = %s\n', k, char(c_vals(k)));
        end
        xlocal = fundamentals * c_vals;
        fprintf('\nParticular solution x(t):\n');
        disp(xlocal);
    end

    if nargout >= 1, x = xlocal; end
end


function chain = build_real_chain(V, bs, bsize, lam, t)
    % Real-eigenvalue chain: x_j(t) = e^{lam t} * sum_{i=0}^{j-1} (t^i/i!) v_{j-i}.
    n = size(V, 1);
    chain = sym(zeros(n, bsize));
    for j = 1:bsize
        expr = sym(zeros(n, 1));
        for i = 0:j-1
            expr = expr + (t^i / factorial(sym(i))) * V(:, bs + j - 1 - i);
        end
        chain(:, j) = simplify(expr * exp(lam * t));
    end
end


function cj_b = find_conjugate_block(block_lambdas, used, lam)
    cj_b = [];
    for k = 1:numel(block_lambdas)
        if used(k), continue; end
        if isAlways(simplify(block_lambdas(k) - conj(lam)) == 0, ...
                    'Unknown', 'false')
            cj_b = k;
            return
        end
    end
end


function v = clear_rational_fractions(v, t)
    % Multiply column v by the LCM of constant rational denominators of
    % its entries, leaving t-dependent denominators alone. Result is
    % equivalent up to a nonzero scalar (fundamental solutions are
    % defined up to scaling), but free of integer fractions.
    [~, D] = numden(v);
    L = sym(1);
    for k = 1:numel(D)
        d = D(k);
        if has(d, t), continue; end
        try
            dval = double(d);
        catch
            continue;
        end
        if ~isfinite(dval) || abs(dval - round(dval)) > 1e-9, continue; end
        dabs = abs(round(dval));
        if dabs <= 1, continue; end
        L = lcm(L, sym(dabs));
    end
    if ~isequal(L, sym(1))
        % Plain multiplication — symbolic auto-eval distributes through the
        % column. Avoid simplify() here: it would fold trig sums like
        % cos(t)+sin(t) into √2·sin(t+π/4) (phase form), which the Chapter 7
        % slides deliberately avoid.
        v = v * L;
    end
end


function [reals, imags] = build_complex_pair(V, bs, bsize, alpha, beta, t)
    % For lam = alpha + i*beta and a chain V(:, bs:bs+bsize-1),
    % build polynomial chains P_j (real part) and Q_j (imag part),
    % then emit
    %   real_j = e^{alpha t} (P_j cos(beta t) - Q_j sin(beta t))
    %   imag_j = e^{alpha t} (Q_j cos(beta t) + P_j sin(beta t)).
    n = size(V, 1);
    reals = sym(zeros(n, bsize));
    imags = sym(zeros(n, bsize));
    for j = 1:bsize
        P = sym(zeros(n, 1));
        Q = sym(zeros(n, 1));
        for i = 0:j-1
            vc = V(:, bs + j - 1 - i);
            coef = t^i / factorial(sym(i));
            P = P + coef * simplify(real(vc));
            Q = Q + coef * simplify(imag(vc));
        end
        % No simplify(): keep trig in textbook "P cos(βt) − Q sin(βt)" form.
        % simplify() would fold rows like cos(t)−sin(t) into √2·cos(t+π/4),
        % which the Chapter 7 slides explicitly avoid. expand() would unfold
        % sin(2t) into 2·sin(t)·cos(t), also undesired. Symbolic auto-eval
        % already handles exp(0·t)=1 and matrix-scalar distribution.
        reals(:, j) = exp(alpha*t) * (P*cos(beta*t) - Q*sin(beta*t));
        imags(:, j) = exp(alpha*t) * (Q*cos(beta*t) + P*sin(beta*t));
    end
end

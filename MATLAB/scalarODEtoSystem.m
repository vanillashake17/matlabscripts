function A = scalarODEtoSystem(coeffs)
% scalarODEtoSystem - companion matrix for a scalar n-th order linear ODE.
%
% Converts the homogeneous scalar ODE
%   c_n y^(n) + c_{n-1} y^(n-1) + ... + c_1 y' + c_0 y = 0
% into the first-order system Y' = A*Y where Y = [y; y'; y''; ...; y^(n-1)].
% The output A can be passed straight to solveLinearODE.
%
% Useful for past-paper questions that explicitly require this conversion,
% e.g. AY2324 final Q16: "Convert y'' - 6y' + 9y = 0 to a first-order
% system, then solve" (other methods earn 0 marks).
%
% SYNTAX:
%   A = scalarODEtoSystem(coeffs)
%
% INPUT:
%   coeffs - vector of coefficients in HIGHEST-ORDER-FIRST order
%            [c_n, c_{n-1}, ..., c_1, c_0],
%            matching how the equation is written on paper.
%
% OUTPUT:
%   A - n-by-n companion matrix (symbolic) so that Y' = A*Y.
%
% EXAMPLE:
%   A = scalarODEtoSystem([1 -6 9]);   % y'' - 6y' + 9y = 0
%   solveLinearODE(A, [1; 4], 0);       % IVP y(0)=1, y'(0)=4 (AY2324 Q16)

    coeffs = sym(coeffs(:).');  % symbolic row vector
    if numel(coeffs) < 2
        error('Need at least 2 coefficients (a first-order ODE has 2).');
    end
    if isAlways(coeffs(1) == 0)
        error('Leading coefficient (of y^(n)) must be nonzero.');
    end

    n = numel(coeffs) - 1;  % order of ODE

    fprintf('--- Scalar ODE -> first-order system (companion form) ---\n');
    fprintf('ODE: %s = 0\n', pretty_ode(coeffs));

    % Normalise so leading coefficient is 1.
    if ~isAlways(coeffs(1) == 1)
        fprintf('Dividing through by %s:\n', char(coeffs(1)));
        coeffs = coeffs / coeffs(1);
        fprintf('     %s = 0\n', pretty_ode(coeffs));
    end

    fprintf('Substitution: Y = [y; y''; y''''; ...; y^(%d)]\n', n-1);
    fprintf('Then  y^(n) = -c_{n-1} y^(n-1) - ... - c_1 y'' - c_0 y,\n');
    fprintf('and   Y'' = A*Y  with companion matrix:\n');

    % Build n-by-n companion matrix.
    A = sym(zeros(n));
    for i = 1:(n-1)
        A(i, i+1) = sym(1);
    end
    % Last row stores -c_0, -c_1, ..., -c_{n-1}.
    % After normalisation coeffs = [1, c_{n-1}, c_{n-2}, ..., c_1, c_0].
    for j = 1:n
        A(n, j) = -coeffs(n + 2 - j);
    end

    disp(A);

    fprintf('Next:  solveLinearODE(A)                       %% general solution\n');
    fprintf('   or  solveLinearODE(A, [y0; y1; ...; y%d_0], t0)  %% IVP\n', n-1);
end

% -------------------------------------------------------------------------
function s = pretty_ode(coeffs)
% Render a scalar ODE in compact human-readable form.
    n = numel(coeffs) - 1;
    parts = {};
    for k = 1:numel(coeffs)
        ord = n - k + 1;
        c = coeffs(k);
        try
            if isAlways(c == 0)
                continue;
            end
        catch
        end
        if ord == 0
            term = 'y';
        elseif ord == 1
            term = "y'";
        elseif ord == 2
            term = "y''";
        else
            term = sprintf("y^(%d)", ord);
        end

        cstr = '';
        if isAlways(c == 1)
            cstr = '';
        elseif isAlways(c == -1)
            cstr = '-';
        else
            cstr = sprintf('(%s)*', char(c));
        end
        parts{end+1} = sprintf('%s%s', cstr, term); %#ok<AGROW>
    end
    if isempty(parts)
        s = '0';
    else
        s = strjoin(parts, ' + ');
    end
end

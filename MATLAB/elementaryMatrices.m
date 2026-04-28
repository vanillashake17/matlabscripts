function [Es, Eprod, EprodInv] = elementaryMatrices(ops, n)
    % Build the elementary matrices for a sequence of row operations on an
    % n x n identity, plus the cumulative product (so that
    %    Eprod * A   == result of applying ops in order to A).
    % Also returns EprodInv  =  E1^-1 * E2^-1 * ... * Ek^-1 (so that
    %    EprodInv * (Eprod * A) == A).
    %
    % ops is a cell array; each entry is itself a cell describing one op:
    %   {'add',   i, j, k}   -- R_i := R_i + k*R_j
    %   {'swap',  i, j}      -- R_i <-> R_j
    %   {'scale', i, k}      -- R_i := k * R_i   (k must be nonzero)
    %
    % Example (Q2 of midterm: R3+2R1, R1<->R2, -2R2, R3-R2):
    %   ops = { {'add',3,1,2}, {'swap',1,2}, {'scale',2,-2}, {'add',3,2,-1} };
    %   [Es, E, Einv] = elementaryMatrices(ops, 3);

    Es       = cell(1, numel(ops));
    Eprod    = sym(eye(n));
    EprodInv = sym(eye(n));
    for t = 1:numel(ops)
        op = ops{t};
        E    = sym(eye(n));
        Einv = sym(eye(n));
        switch op{1}
            case 'add'
                i = op{2}; j = op{3}; k = sym(op{4});
                E(i, j)    = k;
                Einv(i, j) = -k;
            case 'swap'
                i = op{2}; j = op{3};
                E([i, j], :)    = E([j, i], :);
                Einv            = E;
            case 'scale'
                i = op{2}; k = sym(op{3});
                if isAlways(k == 0, 'Unknown', 'false')
                    error('Cannot scale a row by 0.');
                end
                E(i, i)    = k;
                Einv(i, i) = 1/k;
            otherwise
                error('Unknown row operation: %s', op{1});
        end
        Es{t} = E;
        Eprod    = E * Eprod;
        EprodInv = EprodInv * Einv;
        fprintf('Step %d (%s):\n', t, op{1}); disp(E);
    end
    fprintf('Composite  E_k * ... * E_1  (applies all ops to A):\n'); disp(Eprod);
    fprintf('Inverse    E_1^-1 * ... * E_k^-1  (recovers A from result):\n'); disp(EprodInv);
end

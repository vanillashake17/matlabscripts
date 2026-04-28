function pivots = pivot_columns(R)
% pivot_columns - return the pivot column indices of an RREF matrix.
%
% MATLAB's symbolic `rref` does not support the [R, pivots] = rref(...) two-
% output form (it errors "Too many output arguments"), so any script that
% needs the pivot list for a symbolic matrix has to derive it from R itself.
% This helper does that and works for both numeric and symbolic R.
%
% A row's pivot column is the leftmost nonzero entry. Rows with all zeros
% contribute nothing.
%
% SYNTAX:
%   pivots = pivot_columns(R)
%
% INPUT:
%   R - an RREF matrix (numeric or sym).
%
% OUTPUT:
%   pivots - 1-by-r row vector of pivot column indices, where r = rank(R).

    [rows, cols] = size(R);
    pivots = zeros(1, 0);
    col_start = 1;
    for r = 1:rows
        for c = col_start:cols
            entry = R(r, c);
            try
                isz = isAlways(entry == 0);
            catch
                isz = (entry == 0);
            end
            if ~isz
                pivots(end+1) = c; %#ok<AGROW>
                col_start = c + 1;
                break;
            end
        end
    end
end

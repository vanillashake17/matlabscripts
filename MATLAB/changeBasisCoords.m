function coordsT = changeBasisCoords(coordsS, S, T)
% changeBasisCoords - convert a coordinate vector between two bases of the
% same subspace.
%
% Given [w]_S (coordinates of w with respect to basis S), returns [w]_T such
% that  T*[w]_T = w = S*[w]_S.  Works in exact symbolic arithmetic.
%
% Useful for past-paper questions of the form "given [w]_S = (3,1,2), find
% [w]_T" (AY2324 Q17(b), AY2122 Q1(e), AY2021 Q3(c), AY2223 Q19(b)).
%
% SYNTAX:
%   coordsT = changeBasisCoords(coordsS, S, T)
%
% INPUTS:
%   coordsS - column vector [w]_S (length k = number of columns of S).
%   S       - n-by-k basis matrix; columns span the subspace V.
%   T       - n-by-k alternative basis matrix; columns also span V.
%
% OUTPUT:
%   coordsT - [w]_T satisfying T*coordsT = S*coordsS.
%
% Errors out if T is not a basis (rank-deficient) or w fails to lie in
% span(T) (which means S, T are not bases of the same subspace).

    coordsS = sym(coordsS(:));
    S = sym(S);
    T = sym(T);

    if size(S, 2) ~= numel(coordsS)
        error('Length of coordsS (%d) must equal number of columns of S (%d).', ...
            numel(coordsS), size(S, 2));
    end
    if size(S, 1) ~= size(T, 1)
        error('S and T must have the same number of rows.');
    end

    fprintf('--- Change of coordinates [w]_S -> [w]_T ---\n');

    % Step 1: rebuild w in standard coordinates.
    w = simplify(S * coordsS);
    fprintf('w = S * [w]_S =\n');
    disp(w);

    % Step 2: solve T * coordsT = w via symbolic RREF on [T | w].
    nT = size(T, 2);
    aug = [T, w];
    R = rref(aug);
    pivots = pivot_columns(R);

    if any(pivots == nT + 1)
        error('w is not in span(T). T cannot be a basis containing w.');
    end
    if numel(pivots) < nT
        error('T does not have full column rank: it is not a basis. Check inputs.');
    end

    coordsT = sym(zeros(nT, 1));
    for i = 1:numel(pivots)
        col = pivots(i);
        if col <= nT
            coordsT(col) = R(i, end);
        end
    end

    fprintf('Coordinates [w]_T (so that T * [w]_T = w):\n');
    disp(coordsT);

    % Verify.
    residual = simplify(T * coordsT - w);
    if all(isAlways(residual == 0))
        fprintf('Verified: T * [w]_T = w.\n');
    else
        warning('changeBasisCoords:verifyFailed', ...
            'Verification failed; T * [w]_T does not equal w symbolically.');
    end
end


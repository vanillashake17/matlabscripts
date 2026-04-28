function [wp, wn] = decomposeProjection(w, V)
    % Decompose w into wp + wn, where wp is the orthogonal projection of w
    % onto Col(V) and wn is the orthogonal component (wn perpendicular to V).
    %
    % Uses the projection formula  wp = V (V^T V)^{-1} V^T w.
    % Requires the columns of V to be linearly independent.

    % Always promote to symbolic so machine-epsilon noise (e.g. 1/2^52
    % artefacts from double-precision (V'V)\(V'w)) doesn't leak in.
    V = sym(V);
    w = sym(w(:));
    wpLocal = simplify(V * ((V.' * V) \ (V.' * w)));
    wnLocal = simplify(w - wpLocal);

    fprintf('Projection wp onto Col(V):\n');
    disp(wpLocal);
    fprintf('Orthogonal component wn = w - wp (should satisfy V^T wn = 0):\n');
    disp(wnLocal);

    if nargout >= 1, wp = wpLocal; end
    if nargout >= 2, wn = wnLocal; end
end

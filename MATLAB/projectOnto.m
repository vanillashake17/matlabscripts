function [proj_w_onto_u] = projectOnto(w, u)
%PROJECTONTO Orthogonal projection of w onto u, with full workings printed.
%
%   proj_w_onto_u = projectOnto(w, u)
%
%   Computes  proj_u(w) = ( <w,u> / <u,u> ) * u  and prints every
%   intermediate quantity (the two dot products, the scalar, and the
%   resulting vector) in exact symbolic form.

    w_sym = sym(w(:));
    u_sym = sym(u(:));

    fprintf('\n=== Projection of w onto u ===\n');
    fprintf('w = '); disp(w_sym.');
    fprintf('u = '); disp(u_sym.');

    if all(u_sym == 0)
        fprintf('u = 0  =>  projection is the zero vector.\n\n');
        proj_w_onto_u = zeros(size(w));
        return;
    end

    wu  = simplify(dot(w_sym, u_sym));   % <w,u>
    uu  = simplify(dot(u_sym, u_sym));   % <u,u> = ||u||^2
    coeff = simplify(wu / uu);           % scalar coefficient
    proj  = simplify(coeff * u_sym);     % projection vector

    fprintf('Step 1: dot products\n');
    fprintf('  <w, u> = %s\n', char(wu));
    fprintf('  <u, u> = ||u||^2 = %s\n', char(uu));

    fprintf('\nStep 2: scalar coefficient\n');
    fprintf('  c = <w,u> / <u,u> = %s\n', char(coeff));

    fprintf('\nStep 3: projection vector\n');
    fprintf('  proj_u(w) = c * u =\n');
    disp(proj);

    if isa(w, 'sym') || isa(u, 'sym')
        proj_w_onto_u = proj;
    else
        proj_w_onto_u = double(proj);
    end
    fprintf('================================\n\n');
end

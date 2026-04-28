A = [1 2 3; 4 5 6;]
A(1,:) = A(1,:) + A(2,:)% row addition
A([1,2],:)= A([2,1],:) %row swap
[L,U] = lu(A) %auto LU factorisation
subs(R, a, 2) % sub a into R, temporary unless assigned
subs(R, b, 0)
R_sub = subs(R, [a, b], [2, 5]) % multiple
U([1 2 3],:) %slice first 3 rows
adjoint(A)
det(A)
null(sym(L))% find basis for nullspace
A.' % transpose
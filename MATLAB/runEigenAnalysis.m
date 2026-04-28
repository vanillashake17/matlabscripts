% This script demonstrates the use of the general eigenAnalysis.m function.
clc;
clear;

% --- Case 1: A non-symmetric, diagonalizable matrix ---
fprintf('=============== CASE 1: NON-SYMMETRIC, DIAGONALIZABLE ===============\n');
A1 = [4, -2; 1, 1];
disp('Input Matrix A1:');
disp(A1);
[P1, D1] = eigenAnalysis(A1);

% --- Case 2: A non-diagonalizable matrix ---
fprintf('\n\n=============== CASE 2: NON-DIAGONALIZABLE ===============\n');
A2 = [5, 1; 0, 5];
disp('Input Matrix A2:');
disp(A2);
[P2, D2] = eigenAnalysis(A2);

% --- Case 3: A 3x3 diagonalizable matrix ---
fprintf('\n\n=============== CASE 3: 3x3 DIAGONALIZABLE ===============\n');
A3 = [ 2, -1,  1;
       0,  3, -1;
       2,  1,  3];
disp('Input Matrix A3:');
disp(A3);
[P3, D3] = eigenAnalysis(A3);

% --- Case 4: A symmetric matrix (will trigger orthogonal diagonalization) ---
fprintf('\n\n=============== CASE 4: SYMMETRIC MATRIX ===============\n');
A4 = [ 2, -1, 0;
      -1,  2, -1;
       0, -1, 2];
disp('Input Matrix A4:');
disp(A4);
[P4, D4] = eigenAnalysis(A4);


% The resulting variables are now in the workspace.
disp(' ');
disp('Script finished. Check the workspace for the resulting matrices.');
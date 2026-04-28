% This script runs the svd_exact function on the matrix from the
% MA1522 Week 12 Problem Solving Session slides.

clc;
clear;

% Define the matrix from Problem 1(a)
A = [1  0  1  0;
    -1  2  1  2;
     1  1  2  1];

% Run the exact SVD function
[U_matrix, S_matrix, V_matrix] = svd_exact(A);

% The final results U_matrix, S_matrix, and V_matrix are now available
% in the workspace for further use.
disp(' ');
disp('Script finished. Check the workspace for the final matrices:');
disp('U_matrix, S_matrix, V_matrix');
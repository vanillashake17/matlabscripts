% Clear standard outputs
clc; clear;

% Define the input vectors
u1 = [1; 0; 1; -1; 1];
u2 = [1; 0; 1; 0; 0];
u3 = [1; 1; -1; 1; -1];

% Combine them into a single matrix (columns are the vectors)
U_input = [u1, u2, u3];

disp('--------------------------------------------------------');
disp('Starting Rational Orthonormalization with strict checks');
disp('--------------------------------------------------------');

% Call the function
% This will display the rational result and then error out
E_rational = orthonormalize_rational(U_input);

% (This line will only be reached if the input vectors were already orthogonal)
disp('Process complete without errors.');
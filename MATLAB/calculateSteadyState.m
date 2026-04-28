function steadyStateColumnVector = calculateSteadyState(P)
%calculateSteadyState Exact symbolic steady-state vector of a Markov chain.
% Computes pi such that P*pi = pi, ||pi||_1 = 1, using symbolic arithmetic.
%
% Input:
%   P: A square transition probability matrix (can be numeric or symbolic).
%
% Output:
%   steadyStateColumnVector: A symbolic column vector representing the exact
%                            steady-state distribution.

% Step 1: Validate the input and convert to symbolic form
fprintf('--- Step 1: Validating and converting the transition matrix to symbolic form ---\n');
[n, m] = size(P);
if n ~= m
    error('The transition matrix must be a square matrix.');
end
fprintf('Matrix is square with %d states.\n', n);

% Convert the numerical matrix to a symbolic matrix for exact calculations
P_sym = sym(P);
fprintf('Symbolic transition matrix P:\n');
disp(P_sym);

% Step 2: Set up the symbolic system of equations
fprintf('\n--- Step 2: Setting up the symbolic system of equations ---\n');
% We need to solve the system pi*P = pi, which is equivalent to (P' - I)*pi' = 0.
% Let pi be the unknown steady-state column vector.
pi = sym('pi', [n 1]);
assume(pi, 'real'); % Assume the components are real numbers

% Create the identity matrix of the same size
identity = eye(n);

% Form the system of equations (P' - I) * pi = 0
system_lhs = P_sym - identity;
system_equations = system_lhs * pi == zeros(n, 1);

fprintf('The system is defined by (P'' - I) * pi = 0:\n');
disp(system_equations);

% Step 3: Add the constraint that probabilities must sum to 1
fprintf('\n--- Step 3: Adding the constraint that the sum of probabilities is 1 ---\n');
% This system is linearly dependent. We will replace the last equation
% with the constraint that the sum of the elements of pi must equal 1.
sum_constraint = sum(pi) == 1;
fprintf('Sum constraint: sum(pi) = 1\n');

% Combine the first n-1 equations from the original system with the sum constraint
final_system = [system_equations(1:n-1); sum_constraint];
fprintf('Final system of equations to be solved:\n');
disp(final_system);

% Step 4: Solve the system for the steady-state vector
fprintf('\n--- Step 4: Solving the system of equations symbolically ---\n');
% The 'solve' function finds the exact solution.
solution = solve(final_system, pi);

% The solution is returned as a structure. We need to extract the values.
field_names = fieldnames(solution);
steadyStateColumnVector = sym(zeros(n, 1));
for i = 1:n
    steadyStateColumnVector(i) = solution.(field_names{i});
end

fprintf('Exact steady-state vector (in column form):\n');
disp(steadyStateColumnVector);

% Step 5: Verification (optional)
fprintf('\n--- Step 5: Verifying the exact solution ---\n');
% We check if the result satisfies pi' * P = pi' (or P' * pi = pi)
verification = P_sym * steadyStateColumnVector;
fprintf('Result of P'' * pi:\n');
disp(verification);

% Check if the sum is exactly 1
sum_val = sum(steadyStateColumnVector);
fprintf('Sum of vector elements: %s\n', char(sum_val));

end
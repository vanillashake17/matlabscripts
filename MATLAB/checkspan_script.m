fprintf('Running examples...\n\n');

% Example 1: Equal Spans
disp('--- Example 1: Equal Spans ---');
A1 = [1 0; 0 1; 0 0];
B1 = [1 2; 3 4; 0 0];
compare_spans(A1, B1);

% Example 2: Proper Subset (A is a subset of B)
disp('--- Example 2: A is a proper subset of B ---');
A2 = [3; 4; 0];
B2 = [1 0; 0 1; 0 0];
compare_spans(A2, B2);
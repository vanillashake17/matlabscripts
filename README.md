# MA1508E

NUS MA1508E (Linear Algebra for Engineering) — a small MATLAB toolbox of symbolic helpers for the course, plus a single combined reference doc.

Course material (slides, tutorials, homework, midterms, past papers) lives alongside the code locally but is gitignored.

## What's in here

- `MATLAB/` — symbolic-arithmetic `.m` scripts. One script per skill: RREF, Gram–Schmidt, eigen analysis, SVD, ODE solver, projection, least squares, etc. Designed to print human-readable workings to the command window so they double as exam aids.
- `Summaries/MATLAB_Reference.md` — the canonical reference. Quick-lookup table grouped by chapter, plus a worked example for each core function. The matching `.pdf` is generated from this markdown.

## Using the toolbox

1. Open MATLAB and add the `MATLAB/` folder to the path (or open `MATLAB/MATLAB.prj`).
2. Call any function directly, e.g.

   ```matlab
   solveLinearSystem_RREF(A, b)
   [P, D] = eigenAnalysis(A)
   gramSchmidt(V)
   ```

3. All scripts default to **exact symbolic arithmetic** — fractions and surds, no decimals. Pass numeric arrays in; the function wraps them in `sym(...)` internally.

See `Summaries/MATLAB_Reference.md` for the full list.

## Disclaimer

These scripts are study and revision aids. They've been spot-checked against tutorial and past-paper answers but are not guaranteed correct in every edge case — verify the workings before trusting an answer.

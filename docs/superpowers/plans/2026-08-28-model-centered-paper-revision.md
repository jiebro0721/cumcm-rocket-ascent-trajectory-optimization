# Model-Centered Rocket Paper Revision Implementation Plan

> **For agentic workers:** Execute this plan task by task. Preserve unrelated user changes and keep all paper edits inside `论文/main.tex`.

**Goal:** Correct the optimization model's two material numerical defects, re-evaluate Questions 3 and 4 without brute-force multi-start searches, and rewrite the paper so the mathematical model—not DOP853/SLSQP—is the main argument.

**Architecture:** Retain the repository's shared planar ECI variable-mass dynamics and terminal circular-orbit manifold. Treat each flight phase as a state-transition map, use Question 2's feasible solution as the warm start for higher-dimensional controls, and validate every candidate by independent high-accuracy propagation. Correct throttled fuel feasibility through an integral/mass constraint rather than the full-thrust burn-time bound.

**Tech Stack:** Python 3, NumPy, SciPy (`solve_ivp`, nonlinear equations/constrained optimization), pytest-style numerical assertions, XeLaTeX, Git.

---

### Task 1: Establish the numerical trust boundary

**Files:**
- Inspect: `src/common.py`
- Inspect: `src/q2_inscription.py`
- Inspect: `src/q3_fuel_opt.py`
- Inspect: `src/q4_throttle_opt.py`
- Inspect: `results/*.csv`

**Step 1:** Record the current Question 1–4 outputs and terminal residuals.

Run: `python src/q1_baseline.py` and the lightweight evaluation entry points for Questions 2–4.

**Step 2:** Demonstrate the coast-time derivative loss caused by `round(t_coast, 6)` and the invalid throttled burn-time rejection with targeted perturbations.

**Step 3:** Save the observations in the final comparison; do not treat the current Question 3/4 optima as validated evidence.

### Task 2: Repair the model implementation

**Files:**
- Modify: `src/q2_inscription.py`
- Modify: `src/q3_fuel_opt.py`
- Modify: `src/q4_throttle_opt.py`
- Modify only if required for throttle-aware propagation: `src/common.py`

**Step 1:** Replace rounded coast-state cache keys with exact float keys (or remove the cache) so finite-difference perturbations remain visible.

**Step 2:** For Question 4, replace `t_shut <= t_2+t_{b2}` by the physically correct propellant-capacity condition

\[
\int_{t_2}^{t_3}\sigma(t)\,dt\le t_{b2},
\qquad
m(t_3)\ge m_{s2}+m_{pl}.
\]

Use `t_2+t_{b2}/\sigma_{\min}` only as a loose finite search bound.

**Step 3:** Keep the optimization objective pure (fuel used); use terminal residuals as hard constraints and reserve penalties only for failed propagation.

**Step 4:** Add small regression checks for coast-time sensitivity and throttle-feasible duration.

### Task 3: Recompute Questions 3 and 4 efficiently

**Files:**
- Modify: `src/q3_fuel_opt.py`
- Modify: `src/q4_throttle_opt.py`
- Update generated outputs under: `results/`

**Step 1:** Warm-start Question 3 from the verified Question 2 shooting solution; use nested control grids so a coarse feasible control is representable on the refined grid.

**Step 2:** Compare at least two discretizations using identical constraints. Require the refined optimum not to be worse than the embedded coarse feasible solution within numerical tolerance.

**Step 3:** Warm-start Question 4 from the corrected Question 3 solution with `\sigma=1`; allow longer low-throttle arcs through the corrected fuel constraint.

**Step 4:** Independently re-integrate final candidates at tighter tolerances and report physical terminal errors in metres and metres per second. Avoid large blind multi-start grids.

### Task 4: Rewrite the mathematical narrative

**Files:**
- Modify: `论文/main.tex`

**Step 1:** Rewrite the abstract and problem analysis around the hierarchy

\[
\text{given-control propagation}
\to\text{finite-dimensional shooting}
\to\text{direction optimal control}
\to\text{direction--throttle joint control}.
\]

Remove run counts and solver brand names from the abstract.

**Step 2:** Present the ascent as a hybrid state-transition model with phase maps and a mass jump. Reduce DOP853 to one implementation sentence.

**Step 3:** Correct the terminal circular-orbit manifold to the prograde conditions `r=a`, `v_r=0`, `v_t=sqrt(mu/a)` and use `atan2` for flight-path angle.

**Step 4:** Describe Question 2 as a square nonlinear shooting system and a feasible insertion solution, not an overdetermined system or optimum.

**Step 5:** For Question 3, emphasize the fixed-thrust equivalence between minimum fuel and minimum burn duration; describe control transcription, feasibility, mesh consistency, and the limits of local optimality claims.

**Step 6:** For Question 4, use the thrust-integral propellant constraint and only state conclusions supported by the recomputed solution.

**Step 7:** Add inline citations to the existing bibliography, replace unsupported global-optimality claims, and add a concise conclusion section. Do not edit other LaTeX files.

### Task 5: Verify the artifact and prepare review

**Files:**
- Verify: `论文/main.tex`
- Verify: generated `论文/main.pdf`

**Step 1:** Compile with XeLaTeX twice in a temporary/output-safe manner.

Run from `论文`: `xelatex -interaction=nonstopmode -halt-on-error main.tex` twice.

**Step 2:** Inspect compilation warnings, equation/table references, citations, and representative rendered pages.

**Step 3:** Run numerical regression checks and `git diff --check`.

**Step 4:** Summarize the old-vs-new model comparison, commits, and any residual limitations. Push the feature branch and open a pull request only if repository authentication permits it without prompting.

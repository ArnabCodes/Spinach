# Spinach optimal-control method selection

Date: 2026-04-21

## Question
When should Spinach optimal-control work stay with L-BFGS versus switching to full Newton-GRAPE?

## Sources reviewed
- https://spindynamics.org/wiki/index.php?title=Optimal_control_module
- https://spindynamics.org/wiki/index.php?title=Lbfgs.m
- https://spindynamics.org/wiki/index.php?title=Linesearch.m

## Concrete takeaway
Spinach documents three related optimisation modes for GRAPE: plain GRAPE, LBFGS-GRAPE, and Newton-GRAPE. The practical default is L-BFGS via `fminnewton.m`, because Spinach describes it as "a good mix of computational efficiency and fast convergence" while avoiding explicit Hessian formation or inversion. The line-search helper also explicitly assumes cheap gradients, which fits the quasi-Newton path where many candidate steps may be explored without paying full second-derivative cost.

Full Newton-GRAPE is not the default. The optimal-control module says to request it explicitly by setting `control.method='newton'` instead of `lbfgs`. That implies a selection rule: stay with L-BFGS for general-purpose optimisation and switch to Newton only when second-derivative information is available and the extra curvature cost is justified by the problem structure.

## Working rule for future use
- Start with L-BFGS unless there is a specific reason to exploit exact second-derivative structure.
- Treat L-BFGS as the baseline for large or routine Spinach control problems because it approximates Newton directions without explicit Hessian storage.
- Consider `control.method='newton'` only for problems where Hessian-quality curvature information is available and likely to repay its computational overhead.
- Expect line search to be part of the optimisation loop, and remember Spinach's bundled `linesearch.m` is designed for the cheap-gradient setting.

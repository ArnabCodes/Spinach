# Spinach performance tuning notes

_Last updated: 2026-04-19_

## Why this note exists

A compact reference for future Spinach work when runtime matters, based on a quick external doc review plus local example patterns.

## High-confidence takeaways

1. **Start from the conservative examples, then add acceleration deliberately.**
   The Spinach installation guide says the example set is intentionally conservative, and recommends `sys.enable={'greedy','gpu'};` only when the machine has an FP64-capable Nvidia GPU, more than 16 GB of RAM, and more than 4 CPU cores.

2. **`'gpu'` is best treated as a large-state-space accelerator.**
   Appendix F says GPU support is available in `evolution.m`, `propagator.m`, `slowpass.m`, `krylov.m`, and `step.m`, and that speedups of 10x or more are typical once state spaces exceed about 50,000 dimensions. If GPU memory overflows, reduce parallel worker count. GPU-to-worker binding can be overridden with `sys.gpu_bind=[a b c ...];`.

3. **`'greedy'` is not a universal speed flag.**
   Appendix F says `sys.enable={'greedy'};` allows worker processes to use as much CPU as they want, which helps when the state space is dominated by one large subspace. The same page warns that it is often counterproductive when efficient parallelisation already exists, for example powder averaging. The setting also persists in the Matlab parallel pool until the pool is restarted or a job with different settings is run.

4. **Scratch and parallel-pool setup matter.**
   The installation page says Spinach should live on low-latency storage, not cloud-synchronised folders, and that scratch storage must be writable. Appendix F adds that `sys.scratch` should point to fast storage reachable by all workers because Matlab Distributed Computing Server and Spinach both write heavily to disk. The installation guide also recommends disabling automatic parallel pool shutdown.

5. **Parallel control has multiple layers.**
   Appendix F documents `sys.parallel` for cluster profile plus worker count, `sys.parprops` for cluster `AdditionalProperties`, and notes that these are ignored if a parallel pool is already running. Local examples also show experiment-level controls such as `control.parallel='ensemble';` in optimal-control code.

6. **If Matlab's pool wedges after an interrupt, plan for recovery.**
   The Spinach FAQ says Ctrl-C does not stop workers until their current chunks finish. In mild cases, use `smack.m` to force pool shutdown. If Matlab starts throwing `Composite/subsref` invalid-indexing errors, the documented fix is to restart Matlab.

## Local code patterns worth copying carefully

- `examples/nmr_proteins/hsqc_ubiquitin_a.m` uses `sys.enable={'greedy'}; % 'gpu'` for a large liquid-state protein simulation.
- `examples/nmr_solids/mas_powder_trp_fplanck.m` comments out GPU in one section and enables `sys.enable={'gpu'};` in another, which matches the doc warning that acceleration choices are workload-specific.
- `examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m` uses `control.parallel='ensemble';` for GRAPE optimisation, showing that some parallelisation choices live in experiment controls rather than only in `sys.enable`.

## Practical default for future work

When a future Spinach task looks slow, check these in order:

1. Confirm the machine is suitable: Matlab version, RAM per core, SSD-backed scratch, and GPU class.
2. Identify the dominant workload: one huge subspace, powder/orientation averaging, optimal-control ensemble, or repetitive propagator reuse.
3. Only then choose acceleration:
   - try `{'gpu'}` for large state spaces on supported hardware
   - try `{'greedy'}` when one large subspace dominates
   - avoid assuming `'greedy'` helps powder averaging
   - consider cache switches for repetitive pulse sequences
4. If timings are surprising, restart the parallel pool before benchmarking again because prior pool settings may persist.

## Sources

### External Spinach docs
- https://spindynamics.org/wiki/index.php?title=Installation
- https://spindynamics.org/wiki/index.php?title=Appendix_F:_expert-level_options
- https://spindynamics.org/wiki/index.php?title=Frequently_asked_questions

### Local workspace examples
- `/home/administrator/.openclaw/workspace/Spinach/examples/nmr_proteins/hsqc_ubiquitin_a.m`
- `/home/administrator/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_trp_fplanck.m`
- `/home/administrator/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m`

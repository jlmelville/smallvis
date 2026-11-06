# Native allocation benchmarks

This is a developer evidence script, separate from package validation. It compares separately
installed versions in fresh R processes using fixed synthetic inputs. It needs R, each version of
smallvis with its normal dependencies installed, and GNU `/usr/bin/time` for process peak RSS. It
adds no package dependencies. Baseline means the Phase 4 revision, before removal of the unused
conversion; candidate means the Phase 5 implementation. The script checks the compiled registration
and, for the baseline, exercises the old conversion's distinguishing behavior.

From the repository root, a small harness smoke is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript benchmarks/native-allocations.R /tmp/smallvis-candidate-lib candidate \
  gradient 64 1 2 /tmp/smallvis-gradient-smoke.csv
```

For evidence, install each revision with the same compiler and build flags. Run each command at
least three times per variant, alternating variants, and repeat with `THREADS` set to 1 and 2.
Use unique output filenames for each repetition. For example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time \
  -f '%M' -o /tmp/smallvis-gradient-rss.txt \
  Rscript benchmarks/native-allocations.R /tmp/smallvis-candidate-lib candidate \
  gradient 5000 1 10 /tmp/smallvis-gradient.csv
```

The cases are:

- `gradient`: repeated calls with fixed P, W, Z, and coordinates. P and W are constructed before
  timing. Each call computes a dense native gradient.
- `calibration`: repeated full dense Gaussian calibrations on eight-dimensional coordinates.
- `embedding`: one complete public t-SNE call. `ITERATIONS` sets its optimization budget. Timing
  includes input preparation, optimization, and the final epoch calculation.

The CSV records timed elapsed seconds, time divided by the requested iteration count, checksums,
input sizes, threads, R/platform, and loaded paths. For embedding, time per iteration includes
amortized preparation and final scoring. Thread creation is included; only a two-point native
warmup and garbage collection precede timing. The finite-result scan is included and identical
between variants.

GNU time reports maximum resident memory for the entire process in KiB, including package loading,
data setup, warmup, and computation. It is not an allocation counter or a phase-specific peak.
R's allocation profiler does not see ordinary C++ vector allocations. Interpret process RSS
alongside the removed buffer sizes and the complete embedding measurements.

Write raw CSV and RSS outputs outside the checkout, as above. Keep a small result table and the
machine/compiler details in the active implementation plan. Large generated data and binary outputs
should stay outside Git.

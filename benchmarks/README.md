# Native allocation benchmarks

This is a developer evidence script, separate from package validation. It compares separately
installed versions in fresh R processes using fixed synthetic inputs. It needs R, each version of
smallvis with its normal dependencies installed, GNU `nm`, and GNU `/usr/bin/time` for process peak
RSS. It adds no package dependencies. For the current allocation comparison, baseline is the
completed correctness revision `9e5de8d`; candidate borrows the MDS input matrices and fills Student
weights directly into R storage. The script verifies the distinguishing compiled helper signatures
in the loaded DLL. These revisions intentionally have identical numerical behavior. Use unstripped
Linux builds that retain those symbols; the script fails if it cannot establish their identity.

`REVISION` records the source revision used to build each installation. Keep the build logs, source
file hashes, and compiler settings beside the raw results so the supplied label can be checked.
The CSV additionally records the loaded DLL checksum and helper signatures. Historical timings
made with the former conversion-based variant check have different boundaries and are not directly
comparable to these runs.

From the repository root, a small harness smoke is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla benchmarks/native-allocations.R /tmp/smallvis-candidate-lib candidate \
  REVISION student-weights 64 1 2 /tmp/smallvis-student-smoke.csv
```

For evidence, install each revision with the same compiler and build flags. Run each command at
least three times per variant, alternating variants, and repeat with `THREADS` set to 1 and 2.
Use unique output filenames for each repetition. For example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time \
  -f '%M' -o /tmp/smallvis-gradient-rss.txt \
  Rscript --vanilla benchmarks/native-allocations.R /tmp/smallvis-candidate-lib candidate \
  REVISION mds-gradient 5000 1 3 /tmp/smallvis-gradient.csv
```

The cases are:

- `gradient`: repeated calls with fixed P, W, Z, and coordinates. P and W are constructed before
  timing. Each call computes a dense native gradient.
- `calibration`: repeated full dense Gaussian calibrations on eight-dimensional coordinates.
- `embedding`: one complete public t-SNE call. `ITERATIONS` sets its optimization budget. Timing
  includes input preparation, optimization, and the final epoch calculation.
- `mds-gradient`: repeated native MDS gradients with fixed target/current distances and coordinates.
  Column-wise distance setup precedes timing and avoids matrix-sized conversion scratch storage.
- `student-weights`: repeated native Student-weight calculations on fixed two-dimensional
  coordinates. The native unit diagonal is retained.
- `mds-embedding`: one complete public metric MDS call, with the same timing boundary as
  `embedding`.

The embedding cases use fixed steepest-descent steps and disable tolerance-based early stopping;
they check that the requested iteration budget completed. MDS uses a step of `0.01 / N` and t-SNE
uses `0.1`. These short runs measure the actual consumers, rather than convergence quality.

The CSV records timed elapsed seconds, time divided by the requested iteration count, squared
result checksums, final embedding costs, input sizes, threads, R/platform, BLAS settings, and loaded
paths. For embedding, time per iteration includes amortized preparation and final scoring. Thread
creation is included; only a two-point native warmup and garbage collection precede timing. Finite
checks and checksums follow timing and scan columns to avoid another matrix-sized allocation. A
previous result can remain live while the next call returns, as in an ordinary R assignment loop.
The `dimensions` field records the width of synthetic X (eight); Y always has two columns.

GNU time reports maximum resident memory for the entire process in KiB, including package loading,
data setup, warmup, and computation. It is not an allocation counter or a phase-specific peak.
R's allocation profiler does not see ordinary C++ vector allocations. Interpret process RSS
alongside the removed buffer sizes and the complete embedding measurements.
For Student weights, also run `ITERATIONS=1` to distinguish single-call staging overlap from the
retained outputs and garbage collection of a repeated assignment loop. Keep both results. Even
column-wise post-timing verification contributes allocations to the process peak.

Keep BLAS threads fixed, record compiler flags, and avoid other benchmark/check jobs during evidence
runs. Record host load and any available laptop power information at case boundaries. Investigate
material regressions and report noisy results without asserting a speedup.

Write raw CSV and RSS outputs outside the checkout, as above. Keep a small result table and the
machine/compiler details in the active implementation plan. Large generated data and binary outputs
should stay outside Git.

## Barnes-Hut t-SNE optimizer steps

[`bh-optimization.R`](bh-optimization.R) compares the original t-SNE implementation at
`a3f77cd` with the shared traversal. Install each version into a separate library with
`R CMD INSTALL --preclean` and matching release flags. Run each invocation in a fresh R process.
The script requires R, smallvis and its regular dependencies, and GNU `nm` for the compiled
implementation check. Use unstripped Linux builds.

A small harness smoke from the repository root is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla benchmarks/bh-optimization.R /tmp/smallvis-candidate-lib shared 200 1 2 1
```

The script uses seed 20261107, eight-dimensional Gaussian inputs, and the first two input columns
divided by two as initial coordinates. It builds approximate nearest-neighbor affinities with
perplexity 30 before timing. Each repetition starts from the same state and runs the requested
number of delta-bar-delta steps with `theta = 0.5` and learning rate `N / 12`. One full batch
warms the implementation. The timed region includes tree construction, gradient evaluation, and
the R optimizer update. Input preparation, objective scoring, and final reporting are outside
that region.

For evidence, run each size and thread count with library/variant order
`baseline, shared, shared, baseline`. Use 20 steps and three repetitions per process. The script
writes one CSV row per repetition to standard output with fields
`variant,n,threads,steps,nnz,elapsed,cpu,checksum,dll_md5`. Time columns are in seconds and the
embedding checksum is the sum of squared coordinates after the final step. Write raw results and
build logs outside the repository. This is a developer benchmark, separate from package validation.

The November 7, 2026 comparison used R 4.5.2 on Linux x86-64 with an Intel Core i7-8750H,
GCC 15.2.0 with `-O2`, and one thread for BLAS, OpenMP, and MKL. Each cell reports six batches
of 20 steps per variant. Numbers are median elapsed milliseconds per step, with the sample range
in parentheses.

| Points | Native workers | Original | Shared | Change |
| ---: | ---: | ---: | ---: | ---: |
| 2,000 | 1 | 8.33 (8.10–9.70) | 8.28 (7.85–8.75) | −0.6% |
| 2,000 | 2 | 7.38 (6.50–7.80) | 5.58 (5.50–6.05) | −24.4% |
| 10,000 | 1 | 62.58 (58.15–65.55) | 61.93 (58.65–68.55) | −1.0% |
| 10,000 | 2 | 40.05 (38.55–41.20) | 33.93 (32.15–36.55) | −15.3% |

The final-coordinate checksums and affinity counts matched within each workload. The serial
ranges overlap, and the two-worker runs used less elapsed and CPU time with the shared code.
An initial version showed a serial slowdown at 10,000 points. The final implementation binds
the gradient accumulators once per query and writes each worker's normalizer subtotal after
its traversal. The comparison measures these changes together with the shared traversal.

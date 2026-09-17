# MaMaMIA 0.99.0

* Initial submission.

## Performance

* `correctReadCounts()` no longer calls `glmmTMB::predict()` twice on the full
  window table. The fitted ZINB response is now evaluated once for the distinct
  (GC, subgenome) pairs, which avoids the per-call AD rebuild and random-effect
  re-solve. Results agree with the previous implementation to ~1e-13, and the
  step drops from ~21 s to ~3 s on the packaged `triticum` example. The
  remaining cost is the `glmmTMB` fit itself, which scales with `cores`.

## Reproducibility

* `segments()` now pins the seed of the permutation-based p-values used by
  `DNAcopy::segment()` through a new `seed` argument (default `1L`). CBS was
  previously unseeded, so the same object could yield 22 or 23 segments from one
  call to the next. The caller's RNG state is restored afterwards.

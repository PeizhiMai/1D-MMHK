# Clean reproduction test

The clean-test script is `scripts/clean_test.sh`. It builds the pinned portable
native stack in Docker from a clone, executes independent neutral/charged ED
checks and an end-to-end L4 smoke calculation, and then reproduces all 39
archived L10 cases in a separate **solver-free** image. It also authenticates
the complete archive and both historical presentations. Numerical comparisons
and exact test scope are recorded in the adjacent `validation/` reports.

A separate fresh 10-orbital, 50-twist study is also tested; its status and
comparisons are recorded after completion. A smoke test alone is not evidence
of completing that full campaign. Linux/ARM64 on the local laptop is the tested
backend; native macOS, Windows and Linux/AMD64 are not locally certified by
these records. CI is provided for Linux/AMD64 smoke tests after push.

The tests never mount the original Python environment or qbasis installation,
never require the book PDF and never run student code. Published outputs are
numerically reproducible within tolerances, not necessarily binary-identical.
The existing editable seven-slide deck is preserved byte-for-byte.

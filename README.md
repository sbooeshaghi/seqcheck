# seqcheck

`seqcheck` is a planned Rust CLI/library for validating that SGF read content is consistent with the embedded seqspec header.

Planned scope:
- validate read counts and read lengths against seqspec
- validate region-level constraints where possible
- report mismatches with record-level diagnostics
- support machine-readable and human-readable output

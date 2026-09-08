# Synthetic NM/MD regression tests

These tests generate invented DNA and use no downloaded or individual-level
data. They cover [issue #175](https://github.com/bwa-mem2/bwa-mem2/issues/175).

Build on a supported x86-64 system (choose an architecture supported by the CPU):

```sh
make -j4 arch=avx2 CXX=g++
make -C test CXX=g++ nm_md_test
./test/nm_md_test
python3 test/test_nm_md.py --bwa ./bwa-mem2 -v
```

The unit test exercises both CIGAR-generation entry points, both strands,
the ungapped shortcut and banded alignment, score-only calls, N versus N,
lowercase and IUPAC ambiguity, multiple annotation intervals, partially
overlapping runs, and interval boundaries. It checks 96 NM/MD calculations
and 48 score-only calls without any reference files or external tools.

The five Python integration tests require Python 3.7 or later. They build tiny
indexes and compare native SAM output with NM/MD calculated independently from
the original synthetic reference and input reads. Coverage includes paired
reads, insertions, deletions spanning ambiguous bases, clipping, secondary and
supplementary records, and NM values embedded in SA/XA tags.

Optional independent and before/after comparisons:

```sh
python3 test/test_nm_md.py --bwa ./bwa-mem2 \
    --baseline /path/to/unmodified/bwa-mem2 \
    --samtools /path/to/samtools -v
```

`--baseline` checks that the first 11 SAM fields, alignment scores and all
other tags are unchanged except NM/MD and the NM components of SA/XA. The
reference-without-ambiguity control must have identical complete SAM records.
`--samtools` also checks agreement with stock `samtools calmd`. No BAM library
is required by the Python tests. Temporary synthetic files are removed on exit.

## What the fix changes

BWA's packed reference substitutes A/C/G/T at ambiguous positions. The `.amb`
index file retains the original symbols. CIGAR generation now uses these
annotations when emitting NM and MD, after alignment scoring. In accordance
with the [SAM tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf),
only matching A/C/G/T bases count as matches; ambiguous reference bases count
as differences, and MD contains the actual uppercase reference symbols,
including deleted bases. Reverse-strand MD is in forward-reference orientation.

The internal `bwa_gen_cigar` and `bwa_gen_cigar2` interfaces now receive the
reference metadata (`const bntseq_t *`) in place of its length, so they can
access both the length and ambiguity annotations. All in-tree callers are
updated. Existing indexes can be reused; no index format or CLI change is needed.
Seeding, alignment scoring, CIGAR selection and MAPQ calculation retain their
existing behavior. This fix addresses reported NM/MD, including downstream
decisions made using NM; it does not change how BWA searches ambiguous regions.

## Validation of this change

Validated on Linux x86-64 with GCC 9.4.0, `arch=avx2`, and samtools 1.24:

- The tiny N-reference regression fails on unmodified commit `97978f9`.
- All unit and integration tests pass with the fix, including comparisons
  against the unmodified binary and samtools.
- The unit test also passes AddressSanitizer and UndefinedBehaviorSanitizer
  with leak detection enabled for the changed CIGAR-generation source.
- A separate invented-DNA check with the unmodified LevioSAM2 paired workflow
  at `b2bb023b8403a5bc096cba4f7175e3d575241c52` selects candidate A directly:
  A has NM 2 per pair; B now correctly reports NM 4 instead of 0. No intervening
  BAM-tag correction is needed. This demonstrates the downstream effect of
  corrected counts, without establishing biological truth at an ambiguous locus.

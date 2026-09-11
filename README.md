<p align="center">
  <img src="assets/logo.png" alt="cifi-toolkit" width="220">
</p>

# cifi-toolkit

Toolkit for quality control and downstream processing of **CiFi long reads**.

CiFi combines chromosome conformation capture (3C) with PacBio HiFi sequencing. A single CiFi read can contain multiple proximity-ligated genomic segments. `cifi-toolkit` identifies those segments and converts each multi-contact read into Hi-C-like paired contacts for downstream applications such as genome assembly, phasing, contact-map generation, and scaffolding.

- https://cifi.dennislab.org
- https://voles.dennislab.org

---

## Contents

* [Overview](#overview)
* [Installation](#installation)
* [Quick start](#quick-start)
* [Commands](#commands)

  * [`cifi qc`](#cifi-qc)
  * [`cifi digest`](#cifi-digest)
  * [`cifi contacts`](#cifi-contacts)
  * [`cifi filter`](#cifi-filter)
  * [`cifi enzymes`](#cifi-enzymes)
* [How CiFi reads become paired contacts](#how-cifi-reads-become-paired-contacts)
* [Unique segments and `cifi contacts`](#unique-segments-and-cifi-contacts)

  * [Two representations of the same contacts](#two-representations-of-the-same-contacts)
  * [Segment names](#segment-names)
  * [Alignment filtering](#alignment-filtering)
  * [PA5 output and coordinates](#pa5-output-and-coordinates)
  * [Checking the two routes against each other](#checking-the-two-routes-against-each-other)
* [Digest behavior and options](#digest-behavior-and-options)

  * [Minimum number of segments](#minimum-number-of-segments)
  * [Minimum segment length](#minimum-segment-length)
  * [Restriction-site overhang handling](#restriction-site-overhang-handling)
  * [R2 orientation](#r2-orientation)
  * [Pair naming](#pair-naming)
  * [Large datasets](#large-datasets)
* [Digest report](#digest-report)
* [Restriction enzymes](#restriction-enzymes)
* [Custom restriction sites](#custom-restriction-sites)
* [Input formats](#input-formats)
* [Typical workflow](#typical-workflow)
* [Terminology](#terminology)
* [Citation](#citation)
* [License](#license)

---

# Overview

Unlike conventional paired-end Hi-C, a CiFi read can contain **multiple interacting genomic segments in one long PacBio HiFi read**.

`cifi-toolkit` provides five main functions:

| Command         | Purpose                                                                 |
| --------------- | ----------------------------------------------------------------------- |
| `cifi qc`       | Characterize CiFi reads and estimate digestion/contact yield            |
| `cifi digest`   | Digest CiFi reads in silico; paired-end contacts and/or unique segments |
| `cifi contacts` | Reconstruct pairwise contacts from mapped unique segments (YaHS PA5)    |
| `cifi filter`   | Filter aligned paired contacts by mapping status and MAPQ               |
| `cifi enzymes`  | List built-in restriction enzymes and cut positions                     |

A typical analysis looks like:

```text
PacBio CiFi reads
       │
       ▼
    cifi qc
       │
       ▼
  cifi digest
       │
       ├──────────────────────────────┐
       ▼                              ▼
R1 FASTQ + R2 FASTQ           unique segments FASTQ
       │                       (--segments-out)
       ├──────► hifiasm               │
       │        phasing / assembly    ▼
       ▼                        map once (minimap2)
     mapping                          │
       │                              ▼
       ▼                        cifi contacts
   cifi filter                        │
       │                              ▼
       ├──────► contact maps     YaHS PA5 ──► scaffolding
       └──────► scaffolding
```

---

# Installation

Using pip:

```bash
pip install cifi
```

or Bioconda:

```bash
mamba install bioconda::cifi
```

Check the installation:

```bash
cifi --version
cifi --help
```

---

# Quick start

## 1. Check the CiFi library

```bash
cifi qc reads.bam -e HindIII -o qc_out
```

By default, `cifi qc` samples 10,000 reads and reports properties including:

* read-length distribution
* total sequence analyzed
* GC content
* restriction sites per read
* expected segment sizes
* usable-read fraction
* estimated segments per read
* estimated pairwise contacts

The QC step does **not** generate paired FASTQ files.

---

## 2. Generate paired contacts

```bash
cifi digest reads.bam -e HindIII -o sample --gzip
```

This produces:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
sample_stats.json
sample_digestion_report.html
```

The R1 and R2 files contain pairwise contacts generated from the CiFi reads.

---

## 3. Use the contacts downstream

For example, the generated contacts can be supplied to hifiasm for Hi-C-assisted phasing:

```bash
hifiasm \
    --dual-scaf \
    --h1 sample_R1.fastq.gz \
    --h2 sample_R2.fastq.gz \
    -o assembly \
    hifi_reads.fa
```

The paired reads can also be aligned to an assembly for contact-map generation or scaffolding.

---

## 4. Filter mapped pairs

After alignment:

```bash
cifi filter aligned.bam \
    -o filtered.bam \
    -q 30
```

Both mates must be mapped and meet the requested MAPQ threshold for the pair to be retained.

---

# Commands

For detailed options for any command:

```bash
cifi <command> --help
```

---

## `cifi qc`

`cifi qc` provides a quick characterization of a CiFi library before full digestion.

Basic usage:

```bash
cifi qc reads.bam \
    -e HindIII \
    -o qc_out
```

Sample more reads:

```bash
cifi qc reads.bam \
    -e HindIII \
    -n 50000 \
    -o qc_out
```

Analyze all reads:

```bash
cifi qc reads.bam \
    -e HindIII \
    -n 0 \
    -o qc_out
```

Use a custom recognition site:

```bash
cifi qc reads.bam \
    --site GANTC \
    --cut-pos 1 \
    -o qc_out
```

### QC output

The output directory contains an HTML report and machine-readable results, together with tabular and graphical summaries:

```text
qc_out/
├── qc.html
├── qc.json
├── qc.pdf
├── *.tsv
└── *.png
```

The reports summarize the input reads, enzyme-site distribution, expected segment lengths, and estimated pair yield.

---

## `cifi digest`

`cifi digest` converts multi-contact CiFi reads into paired FASTQ contacts.

Basic usage:

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample
```

Compressed output:

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --gzip
```

FASTQ input is also supported:

```bash
cifi digest reads.fastq.gz \
    -e NlaIII \
    -o sample \
    --gzip
```

For each CiFi read, the command:

1. identifies restriction-enzyme sites
2. performs in-silico digestion
3. recovers CiFi segments
4. optionally strips restriction-site remnants
5. removes segments below the minimum emitted length
6. requires a minimum number of retained segments
7. generates every pairwise combination of retained segments
8. optionally reverse-complements R2
9. writes synchronized R1 and R2 FASTQ records

Only segments from the **same original CiFi read** are paired with one another.

With `--segments-out`, each retained segment is additionally written **once**:

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --gzip \
    --segments-out sample.segments.fastq.gz
```

The segments file is gzip-compressed when its name ends in `.gz` and plain
otherwise; `--gzip` only concerns R1/R2. The R1/R2 output, statistics and
report are the same with or without the option. See
[Unique segments and `cifi contacts`](#unique-segments-and-cifi-contacts).

---

## `cifi contacts`

`cifi contacts` takes alignments of the unique segments and rebuilds every
pairwise contact between the confidently mapped segments of each CiFi read,
writing a [YaHS](https://github.com/c-zhou/yahs) PA5 file.

```bash
minimap2 -t 32 -ax map-hifi assembly.fa sample.segments.fastq.gz \
    | samtools sort -n -@ 8 -o sample.segments.ns.bam

cifi contacts sample.segments.ns.bam \
    -o sample.pa5 \
    -q 1
```

This produces:

```text
sample.pa5
sample_contacts_stats.json
sample_contacts_report.html
```

The BAM must be **grouped by read name** (`samtools sort -n`). A header that
declares coordinate order is refused with a message saying so; a header
without a sort order is streamed, and a read whose segments turn out not to be
contiguous stops the run rather than yielding a partial set of contacts.
Processing is streaming: memory holds the segments of one read at a time.

Options:

```text
-q, --mapq      minimum MAPQ for a segment to take part      [default: 1]
-t, --threads   BAM decompression threads                     [default: 4]
--no-json       skip the statistics file
--no-report     skip the HTML report
```

---

## `cifi filter`

After paired contacts have been aligned, `cifi filter` retains confidently mapped pairs.

```bash
cifi filter aligned.bam \
    -o filtered.bam \
    -q 30
```

The default MAPQ threshold is:

```text
30
```

A pair is removed if:

* the records are not paired
* either mate is unmapped
* either mate has MAPQ below the threshold
* the corresponding mate cannot be recovered

Use multiple threads for BAM I/O:

```bash
cifi filter aligned.bam \
    -o filtered.bam \
    -q 30 \
    -t 8
```

An HTML report and JSON statistics are generated by default.

---

## `cifi enzymes`

List the built-in restriction enzymes with:

```bash
cifi enzymes
```

Currently supported:

| Enzyme  | Recognition site | In-silico cut |
| ------- | ---------------- | ------------- |
| DpnII   | `GATC`           | `↓GATC`       |
| MboI    | `GATC`           | `↓GATC`       |
| Sau3AI  | `GATC`           | `↓GATC`       |
| NlaIII  | `CATG`           | `CATG↓`       |
| HindIII | `AAGCTT`         | `A↓AGCTT`     |

The arrow indicates the cut position used for in-silico digestion.

---

# How CiFi reads become paired contacts

A CiFi read can contain several genomic segments that were proximity-ligated during the 3C experiment.

For example:

```text
CiFi read

┌────────────┬─────────┬────────────────┬───────────┐
│ Segment A  │    B    │   Segment C    │ Segment D │
│  1,850 bp  │  720 bp │    2,300 bp    │  1,100 bp │
└────────────┴─────────┴────────────────┴───────────┘
```

`cifi digest` identifies the restriction-enzyme sites and performs an **in-silico digestion**:

```text
A          B          C          D
1850 bp    720 bp     2300 bp    1100 bp
```

It then generates **all pairwise combinations of segments from the same CiFi read**:

```text
A ─ B
A ─ C
A ─ D
B ─ C
B ─ D
C ─ D
```

For a CiFi read containing `n` retained segments:

```text
pairs = n × (n - 1) / 2
```

For example:

| Segments in one CiFi read | Paired contacts |
| ------------------------: | --------------: |
|                         3 |               3 |
|                         4 |               6 |
|                        10 |              45 |
|                        17 |             136 |

The resulting contacts are written as synchronized R1 and R2 FASTQ files.

---

## CiFi paired reads are variable length

The generated R1 and R2 records are **not fixed-length reads** like conventional Illumina Hi-C reads.

Each mate represents one recovered CiFi segment. Because different segments can have different lengths, the two mates can also have different lengths:

```text
Pair 1
R1 = Segment A = 1850 bp
R2 = Segment B =  720 bp

Pair 2
R1 = Segment A = 1850 bp
R2 = Segment C = 2300 bp
```

This is expected.

What is synchronized between the two FASTQ files is the **pairing**, not the sequence length:

```text
sample_R1.fastq.gz       sample_R2.fastq.gz

pair_0001                pair_0001
pair_0002                pair_0002
pair_0003                pair_0003
   ...                       ...
```

Corresponding R1 and R2 records have the **same read name**, remain in the same order, and represent two segments originating from the same CiFi read.

---

# Unique segments and `cifi contacts`

## Two representations of the same contacts

A CiFi read with `n` retained segments carries `n(n-1)/2` pairwise contacts.
`cifi digest` can write those contacts in two forms:

| Output                | Written as                         | Consumer                                  |
| --------------------- | ---------------------------------- | ----------------------------------------- |
| R1/R2 FASTQ           | every pair, both segments repeated | hifiasm (`--h1`/`--h2`), Hi-C style tools |
| unique segments FASTQ | every retained segment once        | an aligner, then `cifi contacts`          |

hifiasm needs actual paired reads, so the R1/R2 form stays as it is. For
scaffolding, however, an aligner only needs to see each segment once: the
contacts can be reconstructed from the mapped coordinates afterwards. With
the pairs form, a read with `n` segments costs `n(n-1)` alignments (two mates
per pair); with unique segments it costs `n`:

| Segments in one read | Contacts | Records mapped, pairs route | Records mapped, segments route |
| -------------------: | -------: | --------------------------: | -----------------------------: |
|                    4 |        6 |                          12 |                              4 |
|                   10 |       45 |                          90 |                             10 |
|                   17 |      136 |                         272 |                             17 |

`cifi contacts` reports both numbers for the input it saw
(`pair_mates_equivalent` against `segments_seen` in the statistics), so the
saving on a real sample can be read off directly.

The intended scaffolding path is therefore:

```text
CiFi reads ── cifi digest ──┬── R1/R2 ──────────────────────────────► hifiasm
                            └── unique segments ── minimap2 ── samtools sort -n
                                                   ── cifi contacts ── PA5 ──► YaHS
```

Segments come out in native read orientation; `--revcomp-r2` is a paired-FASTQ
option and never reaches the segments file. Overhang stripping and
`--min-segment-len` apply to the segment itself, so what is written here is
exactly what the pairs are built from. Reads that fail `--min-segments`
contribute no segments.

## Segment names

Each segment is named after its read with a suffix the toolkit owns:

```text
<original_read_name>__CIFI_SEG__<k>
```

`k` is the 1-based index of the span the segment occupies between cuts in the
read. Dropped spans leave gaps rather than renumbering their neighbours: a read
whose second span was too short yields `__CIFI_SEG__1`, `__CIFI_SEG__3`,
`__CIFI_SEG__4`. The suffix is parsed from the right, so the original name may
contain any number of `/`, `:` or `_` (it must not itself contain
`__CIFI_SEG__`; `cifi digest` refuses such reads). Names survive FASTQ to
SAM/BAM unchanged, which is what `cifi contacts` relies on; a BAM whose names
do not follow the contract is rejected with a clear error rather than being
grouped wrongly.

The R1/R2 pair names keep their existing form, `<read>_<i>_<j-i-1>` with `i`
and `j` counting **retained** segments from 0. The two schemes coexist and
describe the same segments in different terms. PA5 pair names use the segment
scheme, joining the two segment names:

```text
<read>__CIFI_SEG__1__CIFI_SEG__3
```

`cifi.segment_name()` and `cifi.parse_segment_name()` expose the contract to
Python.

## Alignment filtering

For every read, `cifi contacts` keeps the segments whose **primary** alignment
is mapped with MAPQ at or above `-q`, and emits every pair among them. In
particular:

* unmapped segments, secondary and supplementary records are ignored
* a second primary record for a segment already seen (a concatenation of
  alignments, say) is ignored and reported as an anomaly; the first stands
* segments of different reads are never paired; a read with fewer than two
  usable segments emits nothing
* each unordered pair is written once, in segment order; there are no self pairs

The default MAPQ threshold is `1`, following the CiFi paper and the assembly
pipeline that scaffolded with these contacts before. (`cifi filter` defaults to
30, a choice made for contact maps rather than scaffolding.) YaHS applies its
own `-q` on top of the values written to the file.

The statistics file records alignment records seen, distinct segments and
reads, primary mapped / unmapped / secondary / supplementary counts, segments
below the threshold, usable segments, reads with at least two usable segments,
contacts written, the per-read maxima and means, and the mapping-work
comparison against the pairs route.

## PA5 output and coordinates

Each PA5 row has seven tab-separated columns and there is no header line:

```text
pair_name  contig1  pos1  contig2  pos2  mapq1  mapq2
```

YaHS reads it with `yahs -o out assembly.fa sample.pa5` (positions are taken
verbatim, and `MIN(mapq1, mapq2)` is what its `-q` filters on). A `.gz` output
name compresses the file; YaHS reads either.

The position is the **0-based midpoint of the alignment**:
`pos0 + reference_span / 2`, floored, with the reference span summed over the
`M/=/X/D/N` CIGAR operations. This is not a guess: for a name-sorted BAM, YaHS
(`link.c`, `parse_bam_rec` followed by the link write) computes
`s/2 + e/2 + (s&1 && e&1)` from the 0-based start and exclusive end, which is
this value, and its PA5 reader uses the column as given, clamping it to
`length - 1`, so it is treated as 0-based there too. A PA5 written this way
gives YaHS byte-identical links to the ones it would derive from the
equivalent name-sorted BAM. Strand does not enter into it. The lab script that
previously converted segment tables to pairs wrote 1-based alignment starts
into a different (4DN `.pairs`) format, which is why the choice is documented
here rather than inherited.

## Checking the two routes against each other

`tests/e2e/equivalence.py` runs both routes on the same reads and compares the
contact sets exactly (identities, contigs, positions, MAPQs). It needs
`minimap2` and `samtools` and otherwise skips; `pytest tests/e2e` runs it on a
bundled synthetic reference and read set, where the two routes agree exactly.

On a real sample:

```bash
# new route
cifi digest sample.cifi.bam -e HindIII -o sample --gzip \
    --segments-out sample.segments.fastq.gz
minimap2 -t 32 -ax map-hifi assembly.fa sample.segments.fastq.gz \
    | samtools sort -n -@ 8 -o sample.segments.ns.bam
cifi contacts sample.segments.ns.bam -o sample.pa5 -q 1
wc -l sample.pa5

# both routes side by side, with the mapping-work reduction
python tests/e2e/equivalence.py --workdir e2e \
    --reference assembly.fa --reads sample.cifi.bam --threads 32
```

The reduction is the sum over reads of `n(n-1)` (mates the pairs route maps)
divided by the sum of `n` (segments); `cifi contacts` prints it, and the
statistics file carries it as `mapping_work_reduction`.

One difference on real data does not come from the toolkit. minimap2 breaks
ties between equally scoring placements with a hash of the query name, and
the two routes name the same sequence differently, so a small share of
ambiguous (low-MAPQ) segments lands elsewhere in one route than in the other.
The pairs route even disagrees with itself, since it maps each segment once
per pair, under a different name each time; the segments route places each
segment once. With `--minimap2-args=--no-hash-name` the tie-break no longer
depends on the name and the two routes agree exactly.
`scripts/benchmark_contacts.py` times each step of the segments route and
reports peak memory.

---

# Digest behavior and options

## Minimum number of segments

By default, a CiFi read must produce at least:

```text
3 segments
```

after processing.

Change this with:

```text
-m
--min-segments
```

For example:

```bash
cifi digest reads.bam \
    -e HindIII \
    --min-segments 4 \
    -o sample
```

or:

```bash
cifi digest reads.bam \
    -e HindIII \
    -m 4 \
    -o sample
```

---

## Minimum segment length

The default minimum emitted segment length is:

```text
60 bp
```

Change it with:

```text
-l
--min-segment-len
```

For example:

```bash
cifi digest reads.bam \
    -e HindIII \
    --min-segment-len 100 \
    -o sample
```

The threshold applies to the **sequence actually emitted into R1 or R2**, after any requested overhang stripping.

Therefore:

```bash
--min-segment-len 60
```

guarantees that neither emitted mate is shorter than 60 bp.

---

## Restriction-site overhang handling

By default:

```text
--strip-overhang
```

is enabled.

The restriction-site remnant is removed from segments that begin at an in-silico cut before those segments are paired.

Importantly, stripping is performed on the **segment itself**, not specifically on R1 or R2.

For example:

```text
Segment C
    ↓
overhang stripped
    ↓
processed Segment C
```

Every pair containing Segment C then uses the same processed sequence:

```text
A ─ C
B ─ C
C ─ D
```

To retain the restriction-site remnant:

```bash
cifi digest reads.bam \
    -e HindIII \
    --no-strip-overhang \
    -o sample
```

---

## R2 orientation

By default, R2 keeps the orientation in which the segment occurred in the sequenced CiFi read:

```text
--no-revcomp-r2
```

To reverse-complement R2:

```bash
cifi digest reads.bam \
    -e HindIII \
    --revcomp-r2 \
    -o sample
```

Overhang stripping and R2 reverse complementation are **independent options**.

For example:

```bash
cifi digest reads.bam \
    -e HindIII \
    --strip-overhang \
    --revcomp-r2 \
    -o sample
```

first processes the segment boundary and then reverse-complements R2.

---

## Pair naming

The two mates of a generated pair have the **same FASTQ name**.

Conceptually:

```text
R1
@movie/read_123_0_2
ACGT...

R2
@movie/read_123_0_2
TGCA...
```

Pair identity therefore comes directly from the matching FASTQ record names.

R1 and R2 do **not** use different `/1` and `/2` suffixes.

---

## Large datasets

For large datasets:

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --gzip \
    --fast
```

`--fast` uses streaming statistics to reduce memory requirements.

It affects reporting/statistics only; segmentation and generated paired contacts are processed normally.

---

# Digest report

`cifi digest` generates a detailed HTML report by default.

The report includes information about both the original CiFi data and the resulting paired contacts.

### Input reads

Includes:

* number of reads processed
* read-length distribution
* total input bases
* GC content
* restriction-enzyme sites
* sites per read

### Filtering

Reports reads or segments lost because of:

* insufficient enzyme sites
* insufficient retained segments
* segments below the minimum emitted length

### Segment statistics

Includes distributions of:

* segment length
* segments per passing CiFi read
* enzyme sites per read

### Pair statistics

Includes:

* total paired contacts
* pairs generated per CiFi read
* distribution of pair yield across reads

### Base yield

The report also tracks the flow of sequence bases through processing:

```text
input bases
    │
    ├── bases in filtered reads
    ├── short segments removed
    ├── restriction-site remnants trimmed
    │
    └── bases represented in emitted R1/R2 contacts
```

---

# Restriction enzymes

Built-in enzymes:

| Enzyme  | Recognition site | In-silico cut |
| ------- | ---------------- | ------------- |
| DpnII   | `GATC`           | `↓GATC`       |
| MboI    | `GATC`           | `↓GATC`       |
| Sau3AI  | `GATC`           | `↓GATC`       |
| NlaIII  | `CATG`           | `CATG↓`       |
| HindIII | `AAGCTT`         | `A↓AGCTT`     |

---

# Custom restriction sites

Custom recognition sites can be supplied with:

```text
--site
--cut-pos
```

For example:

```bash
cifi digest reads.bam \
    --site GANTC \
    --cut-pos 1 \
    -o sample
```

Degenerate IUPAC bases are supported:

```text
N R Y W S M K B D H V
```

---

# Input formats

`cifi qc` and `cifi digest` support:

```text
FASTQ
FASTQ.gz
BAM
SAM
CRAM
```

For aligned BAM filtering, use:

```text
cifi filter
```

---

# Typical workflow

```text
                         PacBio CiFi reads
                                │
                                ▼
                            cifi qc
                                │
                       library assessment
                                │
                                ▼
                          cifi digest
                                │
                     in-silico digestion
                                │
                         CiFi segments
                                │
                ┌───────────────┴────────────────┐
                ▼                                ▼
      all pairwise contacts             unique segments FASTQ
                │                          (--segments-out)
    ┌───────────┴───────────┐                    │
    ▼                       ▼                    ▼
R1 FASTQ                R2 FASTQ           minimap2, once
    └───────────┬───────────┘                    │
                │                                ▼
         paired contacts                  samtools sort -n
                │                                │
    ┌───────────┴───────────┐                    ▼
    ▼                       ▼              cifi contacts
 hifiasm                 mapping                 │
phasing / assembly          │                    ▼
                            ▼                YaHS PA5
                       cifi filter               │
                            │                    ▼
                  ┌─────────┴─────────┐     scaffolding
                  ▼                   ▼
             contact maps         scaffolding
```

---

# Terminology

### CiFi read

The original long PacBio HiFi sequence produced from a proximity-ligated CiFi molecule.

### Segment

A sequence recovered from a CiFi read by in-silico digestion at the specified restriction-enzyme sites.

### Pair

Two segments from the same CiFi read emitted together as corresponding R1 and R2 FASTQ records.

### Contact

The pairwise relationship represented by those two segments.

A CiFi read containing multiple segments therefore represents multiple pairwise contacts.

### Unique segment

A retained segment written once by `cifi digest --segments-out`, named
`<read>__CIFI_SEG__<k>`, from whose alignment `cifi contacts` rebuilds the
read's contacts.

---

# Citation

If you use `cifi-toolkit` please cite:

*Single-library chromosome-scale diploid assemblies of vole genomes resolve a species-specific duplication implicated in pair bonding.*
**Cell Genomics. 2026; 101336.**
https://doi.org/10.1016/j.xgen.2026.101336

For the CiFi method, please also cite:

*CiFi: accurate long-read chromosome conformation capture with low-input requirements.*
**Nature Communications. 2025.**
https://doi.org/10.1038/s41467-025-66918-y

---

# License

MIT License

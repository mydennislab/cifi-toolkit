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
  * [Output formats: BED and PA5](#output-formats-bed-and-pa5)
  * [Checking the two routes against each other](#checking-the-two-routes-against-each-other)
* [Molecule tables](#molecule-tables)

  * [Two indices](#two-indices)
  * [The segments table](#the-segments-table)
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

`cifi-toolkit` provides these functions:

| Command           | Purpose                                                                 |
| ----------------- | ----------------------------------------------------------------------- |
| `cifi qc`         | Characterize CiFi reads and estimate digestion/contact yield            |
| `cifi digest`     | Digest CiFi reads in silico; paired-end contacts, unique segments and the segments table |
| `cifi contacts`   | Reconstruct pairwise contacts from mapped unique segments (YaHS BED/PA5) |
| `cifi filter`     | Filter aligned paired contacts by mapping status and MAPQ               |
| `cifi enzymes`    | List built-in restriction enzymes and cut positions                     |

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
       ├──────► contact maps     YaHS BED ──► scaffolding
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

A further optional output describes the molecules rather than the
contacts (see [Molecule tables](#molecule-tables)):

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --segments-out sample.segments.fastq.gz \
    --segments-table sample.segments.tsv.gz
```

* `--segments-table PATH`: one row per retained segment (bgzip TSV) with the
  cut span it came from, its rank among the kept segments, the read's
  length and its coordinates in the read. Works with or without
  `--segments-out`; when both are given the rows and the FASTQ records
  correspond one to one.

The option changes neither R1/R2, the segments FASTQ, the report nor the
existing keys of the statistics file; it adds its own keys under
`parameters`, `results` and `output` when used, plus `command` (the command
line) and `formats` (the name and layout version of the table written) at
the top level.

---

## `cifi contacts`

`cifi contacts` takes alignments of the unique segments and rebuilds every
pairwise contact between the confidently mapped segments of each CiFi read,
writing a contacts file for [YaHS](https://github.com/c-zhou/yahs) in either
of the two text formats it reads:

* **BED** (`--format bed`, or an output name ending in `.bed`/`.bed.gz`): two
  consecutive lines per contact, each carrying the actual aligned span of one
  segment. YaHS normalises coverage from those spans and needs no read
  length. This is the format to scaffold CiFi data with.
* **PA5** (`--format pa5`, or an output name ending in `.pa5`/`.pa5.gz`):
  one line per contact with the alignment midpoints. YaHS rebuilds an
  interval around each midpoint from a single global `--read-length`, which
  fits fixed-length Hi-C reads but not CiFi segments.

```bash
cifi digest sample.cifi.bam -e HindIII -o sample --gzip \
    --segments-out sample.segments.fastq.gz

minimap2 -t 32 -ax map-hifi assembly.fa sample.segments.fastq.gz \
    | samtools sort -n -@ 8 -o sample.segments.ns.bam

cifi contacts sample.segments.ns.bam -o sample.bed -q 1 --format bed
yahs -o scaffolds assembly.fa sample.bed
```

This produces:

```text
sample.bed
sample_contacts_stats.json
sample_contacts_report.html
```

The same command with `-o sample.pa5` (or `--format pa5`) writes PA5, exactly
as before; `--format` overrides the extension, with a warning, since YaHS
picks its parser by extension unless given `--file-type`. An output name with
neither extension needs `--format`: rather than guess, the command stops with
a usage error, as YaHS itself would on that name. The output takes its final
name only once it is complete: a run that fails leaves no partial `.bed` or
`.pa5` behind, and an earlier file of that name is left untouched.

The BAM must be **grouped by read name** (`samtools sort -n`). A header that
declares coordinate order is refused with a message saying so; a header
without a sort order is streamed, and a read whose segments turn out not to be
contiguous stops the run rather than yielding a partial set of contacts.
Processing is streaming: memory holds the segments of one read at a time.

Options:

```text
--format        bed or pa5; taken from the output name when omitted,
                required for any other name
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
                                                   ── cifi contacts --format bed ── BED ──► YaHS
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
describe the same segments in different terms. Contact pair names (PA5 column
1, BED column 4) use the segment scheme, joining the two segment names:

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
contacts written, the per-read maxima and means, the mapping-work comparison
against the pairs route, and the output format. The statistics are the same
whichever format was written; only the file differs.

## Output formats: BED and PA5

Both formats carry the same contacts, contig names and MAPQs, filtered the
same way; neither has a header line, and a `.gz` output name compresses the
file (YaHS reads either). They differ in what YaHS can do with the
coordinates.

### BED

Two consecutive lines per contact, five tab-separated columns each:

```text
contig1  start1  end1  pair_name  mapq1
contig2  start2  end2  pair_name  mapq2
```

`start` is the 0-based reference start of the segment's primary alignment
and `end` its exclusive end from the CIGAR (`M/=/X/D/N` consume reference).
Every span is the segment's own: a 60 bp and a 12 kb segment of the same
read are written with 60 bp and 12 kb intervals. Both lines carry the same
pair name.

YaHS reads it with `yahs -o out assembly.fa sample.bed`. What it does with
the file is fixed by its source (`link.c` and `asset.c`, identical between
the `main` branch and v1.2.2):

* the two lines of a contact must be **consecutive**: the reader holds one
  record and pairs it with the next only if the names match; otherwise the
  held record is dropped and the new one takes its place
  (`link.c`, `dump_links_from_bed_file`, lines 1581-1592 and 1652-1660);
* **identical names** on the two lines are accepted (`asset.c`,
  `is_read_pair`, lines 175-188: equal strings, or strings differing only in
  the character after a trailing `/`), so no `/1` `/2` suffix is needed;
* the columns are read as **contig, start, end, name, MAPQ**
  (`sscanf(line, "%s %u %u %s %hhu", ...)`, `link.c` lines 1585 and 1591);
* the link position is the **midpoint** of each interval,
  `s/2 + e/2 + (s&1 && e&1)`, i.e. `floor((start + end) / 2)`
  (`link.c` lines 1633-1634), which is exactly the value the PA5 row carries;
* the **real `start`/`end` interval** goes into the coverage track used for
  normalisation, clamped to the contig (`link.c` lines 1614-1621, then
  `calc_cov_norms` at 1684).

No read length enters into it. This is why BED is the format meant for
scaffolding CiFi data.

### PA5

One line per contact, seven tab-separated columns:

```text
pair_name  contig1  pos1  contig2  pos2  mapq1  mapq2
```

YaHS reads it with `yahs -o out assembly.fa sample.pa5` (positions are taken
verbatim, and `MIN(mapq1, mapq2)` is what its `-q` filters on).

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

What PA5 cannot carry is the extent of the alignment. YaHS's PA5 reader
halves its `--read-length` (default 150) and takes
`[pos - read_length/2, pos + read_length/2]` as every segment's interval for
coverage normalisation (`link.c`, `dump_links_from_pa5_file`, lines 1726 and
1762-1765). One global length is right for Illumina Hi-C pairs and wrong for
CiFi segments, whose lengths span orders of magnitude. PA5 remains supported
as a general contact table; for YaHS scaffolding, write BED.

Reducing a BED file to midpoints reproduces the PA5 file exactly (same
contacts, coordinates and MAPQs), and YaHS derives byte-identical link
records from the two; only its coverage normalisation differs. The test
suite checks the first, and, when `yahs` is on `PATH`, runs both formats
through it.

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
cifi contacts sample.segments.ns.bam -o sample.bed -q 1   # for yahs

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

# Molecule tables

A CiFi read is one HiFi concatemer of restriction-defined segments: an
ordered, multi-way contact molecule. The pairs of `cifi digest` and the
contacts of `cifi contacts` reduce it to pairwise relations; the segments
table described here keeps the molecule whole, as a data file: which
segments it had in which order, which cut spans were dropped and where each
segment sat in the read. It is a bgzip-compressed TSV file (readable with
any gzip reader, and block-addressable by htslib).

The header layout is `##key=value` metadata lines, one item each, in a
fixed order, then the `#columns:` line naming the columns, then the data
rows. The vocabulary is fixed and given below; a reader splits each `##`
line on its first `=`. The first four keys are `cifi_format`, naming the
format, `format_version`, its layout version (`1`), `tool_version`, the
cifi release, and `command`, the command line that wrote the file. The
remaining keys state the conventions of the columns (coordinate systems,
index meanings, the codes of a column). Nothing in the header varies
between two runs of the same command: no timestamp, no host.

`cifi digest`'s statistics file repeats the provenance: it carries
`formats`, keyed by output (`segments_table`), and gains `command` only when
the table was requested; the file is otherwise unchanged.

## Two indices

Two indices describe a segment's place in its molecule, and they are never
overloaded:

| Index            | Meaning                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------ |
| `span_index`     | The 1-based index of the cut span the segment came from. Dropped spans leave gaps. This is the `k` of the segment name `<read>__CIFI_SEG__<k>`: digest provenance. |
| `retained_index` | The dense 1..N order of the segment among the emitted segments of its molecule: within-molecule order for algorithms. |

A read cut into four spans whose second span was too short has segments with
`span_index` 1, 3, 4 and `retained_index` 1, 2, 3.

## The segments table

`cifi digest --segments-table PATH` writes one row per retained segment, in
read order, in the same order as `--segments-out` writes the FASTQ records.

Header:

```text
##cifi_format=segments
##format_version=1
##tool_version=<cifi version>
##command=<command line>
##read_coordinates=0-based-half-open
##span_index=1-based-original-digest-span
##retained_index=1-based-dense-retained-order
##terminal=T5,T3,I,T5T3
#columns: <the columns below>
```

Columns:

```text
molecule_id  span_index  retained_index  segments_kept  spans_total  read_length
read_start  read_end  span_start  span_end  original_len  processed_len  trimmed_5p
terminal  enzyme  cut_offset
```

* `molecule_id`: the read name (no `__CIFI_SEG__` suffix).
* `span_index`, `retained_index`: as above. `segments_kept`: the number of
  retained segments of the molecule (the largest `retained_index`).
  `spans_total`: the number of cut spans of the read, sites + 1, empty spans
  included (the largest possible `span_index`).
* `read_length`: the length of the whole read, repeated on each of its
  rows; with `read_start` and `read_end` it places the segment in the read
  the way a PAF of the molecule would.
* `read_start`, `read_end`: 0-based, half-open coordinates of the emitted
  (post-trim) segment in the original read; `read[read_start:read_end]` is
  the FASTQ record's sequence. `span_start`, `span_end`: the cut span before
  the 5' site-remnant trim. `original_len = span_end - span_start`,
  `processed_len = read_end - read_start`, `trimmed_5p = read_start -
  span_start` (0 for the read's leading span, which begins at no cut, and
  with `--no-strip-overhang`).
* `terminal`: `T5` for `span_index` 1, `T3` for `span_index == spans_total`,
  `I` otherwise. A read of a single span (possible with `-m 1`) is `T5T3`.
* `enzyme`: the enzyme name, or the site itself for `--site`; `cut_offset`:
  the cut position used.

`--gzip` does not change the table; it is always bgzip. Reads that failed
`--min-segments` have no rows.

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
                            ▼            YaHS BED (or PA5)
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

### Molecule

The CiFi read seen as an ordered set of segments: what the segments table
describes.

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

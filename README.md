<p align="center">
  <img src="assets/logo.png" alt="cifi-toolkit" width="220">
</p>

# cifi-toolkit

Toolkit for quality control and downstream processing of **CiFi long reads**.

CiFi combines chromosome conformation capture (3C) with PacBio HiFi sequencing. A single CiFi read can contain multiple proximity-ligated genomic segments. `cifi-toolkit` identifies those segments and converts each multi-contact read into pairwise contacts for downstream applications such as genome assembly, phasing, contact-map generation, and scaffolding.

* https://cifi.dennislab.org
* https://voles.dennislab.org

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
* [Unique segments and native contact reconstruction](#unique-segments-and-native-contact-reconstruction)

  * [Two representations of the same CiFi contacts](#two-representations-of-the-same-cifi-contacts)
  * [Segment names](#segment-names)
  * [Mapping unique segments](#mapping-unique-segments)
  * [Alignment filtering](#alignment-filtering)
  * [BED output](#bed-output)
  * [PA5 output](#pa5-output)
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

| Command         | Purpose                                                                      |
| --------------- | ---------------------------------------------------------------------------- |
| `cifi qc`       | Characterize CiFi reads and estimate digestion/contact yield                 |
| `cifi digest`   | Digest CiFi reads in silico; generate paired contacts and/or unique segments |
| `cifi contacts` | Reconstruct all pairwise contacts from mapped unique segments                |
| `cifi filter`   | Filter aligned paired contacts by mapping status and MAPQ                    |
| `cifi enzymes`  | List built-in restriction enzymes and cut positions                          |

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
       │                              │
       ▼                              ▼
    hifiasm                      minimap2
 phasing/assembly              map each segment once
                                      │
                                      ▼
                               samtools sort -n
                                      │
                                      ▼
                                cifi contacts
                                      │
                                      ▼
                               BED contacts
                                      │
                                      ▼
                                    YaHS
                                  scaffolding
```

The paired R1/R2 representation is useful for tools such as hifiasm that require paired reads.

For scaffolding, the unique-segment route avoids repeatedly mapping the same CiFi segment. Each segment is mapped once and the full set of pairwise contacts is reconstructed afterward from its alignment.

---

# Installation

## Released version

Using pip:

```bash
pip install cifi
```

or Bioconda:

```bash
mamba install bioconda::cifi
```

## Development branch

To install the development branch containing native unique-segment contact reconstruction:

```bash
git clone https://github.com/mydennislab/cifi-toolkit.git
cd cifi-toolkit
git checkout cifi-native-contacts
pip install .
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
cifi qc reads.bam \
    -e HindIII \
    -o qc_out
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

## 2. Generate paired contacts and unique segments

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --gzip \
    --segments-out sample.segments.fastq.gz
```

This produces:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
sample.segments.fastq.gz
sample_stats.json
sample_digestion_report.html
```

The R1 and R2 files contain all pairwise contacts generated from each CiFi read.

The segments file contains each retained CiFi segment **once**.

---

## 3. Use R1/R2 for hifiasm

For example:

```bash
hifiasm \
    --dual-scaf \
    --h1 sample_R1.fastq.gz \
    --h2 sample_R2.fastq.gz \
    -o assembly \
    hifi_reads.fa
```

---

## 4. Map unique segments once

For scaffolding, map the unique segments rather than the pair-expanded R1/R2 files:

```bash
minimap2 -t 32 \
    -ax map-hifi \
    --secondary=no \
    --no-hash-name \
    assembly.fa \
    sample.segments.fastq.gz \
  | samtools sort -n -@ 8 \
      -o sample.segments.ns.bam
```

`samtools sort -n` is required because `cifi contacts` processes all segments belonging to the same original CiFi read together.

`--secondary=no` gives each retained segment one primary placement for contact reconstruction.

`--no-hash-name` prevents query-name-dependent tie breaking for equally scoring mappings.

---

## 5. Reconstruct CiFi contacts

The recommended output for CiFi scaffolding is BED:

```bash
cifi contacts sample.segments.ns.bam \
    --format bed \
    -q 1 \
    -o sample.contacts.bed
```

Then scaffold with YaHS:

```bash
yahs \
    -q 0 \
    --no-contig-ec \
    -o sample.yahs \
    assembly.fa \
    sample.contacts.bed
```

The MAPQ filtering is intentionally applied once:

```text
cifi contacts -q 1
        ↓
      YaHS -q 0
```

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

`cifi digest` performs in-silico restriction digestion of CiFi reads.

Basic paired-contact output:

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

### Unique segment output

With `--segments-out`, each retained segment is additionally written **once**:

```bash
cifi digest reads.bam \
    -e HindIII \
    -o sample \
    --gzip \
    --segments-out sample.segments.fastq.gz
```

The unique segments are the same processed sequences used to generate R1/R2 pairs.

The segments file is gzip-compressed when the output name ends in `.gz`; otherwise it is written as plain FASTQ.

`--gzip` controls R1/R2 compression independently.

The existing R1/R2 output, statistics, and report are unchanged by requesting `--segments-out`.

---

## `cifi contacts`

`cifi contacts` reconstructs pairwise contacts from mapped unique CiFi segments.

Input is a **name-sorted BAM** containing alignments of the segments produced by:

```bash
cifi digest ... --segments-out
```

Recommended mapping:

```bash
minimap2 -t 32 \
    -ax map-hifi \
    --secondary=no \
    --no-hash-name \
    assembly.fa \
    sample.segments.fastq.gz \
  | samtools sort -n -@ 8 \
      -o sample.segments.ns.bam
```

Generate BED contacts:

```bash
cifi contacts sample.segments.ns.bam \
    --format bed \
    -q 1 \
    -o sample.contacts.bed
```

Or PA5:

```bash
cifi contacts sample.segments.ns.bam \
    --format pa5 \
    -q 1 \
    -o sample.contacts.pa5
```

For every original CiFi read containing `k` usable mapped segments, `cifi contacts` emits:

```text
k × (k - 1) / 2
```

contacts.

Processing is streaming: only the segments from one CiFi molecule need to be held in memory at a time.

### Alignment requirements

The BAM must be grouped by read name:

```bash
samtools sort -n
```

Only primary mapped alignments that pass the requested MAPQ threshold participate in contact generation.

Unmapped, secondary, and supplementary alignments are ignored.

Segments from different original CiFi molecules are never paired.

### Common options

```text
--format       output representation: bed or pa5
-q, --mapq     minimum segment MAPQ                [default: 1]
-t, --threads  BAM decompression threads           [default: 4]
--no-json      skip statistics JSON
--no-report    skip HTML report
```

---

## `cifi filter`

After paired R1/R2 contacts have been aligned, `cifi filter` retains confidently mapped pairs.

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

List built-in restriction enzymes with:

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

`cifi digest` identifies restriction-enzyme sites and performs an **in-silico digestion**:

```text
A          B          C          D
1850 bp    720 bp     2300 bp    1100 bp
```

It then generates **all pairwise combinations of retained segments from the same CiFi read**:

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

The contacts can be represented either as paired R1/R2 FASTQ records or reconstructed after mapping the unique segments.

---

## CiFi paired reads are variable length

The generated R1 and R2 records are **not fixed-length reads** like conventional Illumina Hi-C reads.

Each mate represents one recovered CiFi segment.

Because different segments have different lengths, the two mates can also have different lengths:

```text
Pair 1
R1 = Segment A = 1850 bp
R2 = Segment B =  720 bp

Pair 2
R1 = Segment A = 1850 bp
R2 = Segment C = 2300 bp
```

This is expected.

What is synchronized between R1 and R2 is the **pairing**, not sequence length.

Corresponding R1 and R2 records have the same read name, remain in the same order, and represent two segments originating from the same CiFi read.

---

# Unique segments and native contact reconstruction

## Two representations of the same CiFi contacts

A CiFi read with `n` retained segments contains:

```text
n(n-1)/2
```

pairwise contacts.

`cifi digest` can represent those data in two ways:

| Output                | Representation                                                | Typical consumer                          |
| --------------------- | ------------------------------------------------------------- | ----------------------------------------- |
| R1/R2 FASTQ           | every pair explicitly written; segments repeated across pairs | hifiasm `--h1/--h2`, paired-contact tools |
| unique segments FASTQ | every retained segment written once                           | minimap2 → `cifi contacts` → YaHS         |

hifiasm requires paired reads, so R1/R2 remain useful.

For scaffolding, however, repeatedly mapping the same segment is unnecessary.

A read with `n` segments requires:

```text
n(n-1)
```

mapped FASTQ records through the pair-expanded route, because every contact contains two mates.

The unique-segment route requires only:

```text
n
```

mapped records.

For example:

| Segments | Contacts | Pair-route records mapped | Unique segments mapped |
| -------: | -------: | ------------------------: | ---------------------: |
|        4 |        6 |                        12 |                      4 |
|       10 |       45 |                        90 |                     10 |
|       17 |      136 |                       272 |                     17 |

The intended scaffolding path is therefore:

```text
CiFi
 │
 └─ cifi digest
      ├─ R1/R2 ───────────────────────────────► hifiasm
      │
      └─ unique segments
              │
              ▼
         minimap2 map-hifi
              │
              ▼
        samtools sort -n
              │
              ▼
         cifi contacts
              │
              ▼
              BED
              │
              ▼
             YaHS
```

---

## Segment names

Each unique segment is named after its original CiFi read with a toolkit-owned suffix:

```text
<original_read_name>__CIFI_SEG__<k>
```

`k` is the 1-based index of the restriction-defined span occupied by that segment.

Dropped spans leave gaps rather than causing the retained segments to be renumbered.

For example, if the second span is removed because it is too short:

```text
read123__CIFI_SEG__1
read123__CIFI_SEG__3
read123__CIFI_SEG__4
```

This preserves the segment's original position within the digested CiFi molecule.

The suffix is parsed from the right, so normal `/`, `:`, and `_` characters in the original read name are supported.

The original read name must not itself contain:

```text
__CIFI_SEG__
```

Segment names survive FASTQ → SAM/BAM unchanged and allow `cifi contacts` to reconstruct which alignments originated from the same CiFi molecule.

---

## Mapping unique segments

Recommended mapping:

```bash
minimap2 -t 32 \
    -ax map-hifi \
    --secondary=no \
    --no-hash-name \
    assembly.fa \
    sample.segments.fastq.gz \
  | samtools sort -n -@ 8 \
      -o sample.segments.ns.bam
```

### Why `map-hifi`?

CiFi segments are derived from PacBio HiFi reads, so the HiFi preset is appropriate for their base-accuracy profile.

### Why `--secondary=no`?

The standard scaffolding route uses one primary alignment for each segment.

### Why `--no-hash-name`?

minimap2 may otherwise use the query name when breaking ties between equally scoring alignments.

Because pair-expanded and unique-segment representations use different query names for the same sequence, name-dependent tie breaking can lead to different placements for ambiguous segments.

`--no-hash-name` removes that unnecessary representation-dependent behavior.

---

## Alignment filtering

For every CiFi molecule, `cifi contacts` keeps segments whose **primary** alignment is mapped with MAPQ at or above `-q`.

It then emits every unordered pair among those usable segments.

Specifically:

* unmapped records are ignored
* secondary records are ignored
* supplementary records are ignored
* segments below the MAPQ threshold are ignored
* segments from different CiFi molecules are never paired
* a molecule with fewer than two usable segments produces no contacts
* every unordered contact is emitted exactly once

The default `cifi contacts` MAPQ threshold is:

```text
1
```

For the recommended YaHS workflow:

```text
cifi contacts -q 1
        ↓
      YaHS -q 0
```

so MAPQ filtering occurs once.

---

## BED output

BED is the **recommended contact representation for CiFi scaffolding with YaHS**.

Generate it with:

```bash
cifi contacts sample.segments.ns.bam \
    --format bed \
    -q 1 \
    -o sample.contacts.bed
```

Each contact is represented by **two consecutive BED records**, one for each mapped CiFi segment.

Conceptually:

```text
contigA    startA    endA    pair_name    mapqA
contigB    startB    endB    pair_name    mapqB
```

The two records share the same `pair_name`.

Unlike a fixed-read-length representation, BED preserves the actual genomic alignment span of each CiFi segment:

```text
segment A = 350 bp
segment B = 1,800 bp
segment C = 5,200 bp
```

This is useful for CiFi because restriction-defined segments are naturally variable in length.

The BED can be supplied directly to YaHS:

```bash
yahs \
    -q 0 \
    --no-contig-ec \
    -o sample.yahs \
    assembly.fa \
    sample.contacts.bed
```

---

## PA5 output

PA5 is also available when a point-like paired-contact representation is desired:

```bash
cifi contacts sample.segments.ns.bam \
    --format pa5 \
    -q 1 \
    -o sample.contacts.pa5
```

Each PA5 row contains:

```text
pair_name  contig1  pos1  contig2  pos2  mapq1  mapq2
```

The position is the **0-based midpoint of the alignment**.

PA5 represents each segment as a single genomic position rather than preserving its complete alignment span.

For CiFi scaffolding, BED is therefore preferred when using YaHS.

---

# Digest behavior and options

## Minimum number of segments

By default, a CiFi read must produce at least:

```text
3 segments
```

after processing.

Change this using:

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

The threshold applies to the **sequence actually emitted into R1/R2 and the unique-segment output**, after requested overhang stripping.

Therefore:

```bash
--min-segment-len 60
```

guarantees that emitted segments are at least 60 bp long.

---

## Restriction-site overhang handling

By default:

```text
--strip-overhang
```

is enabled.

The restriction-site remnant is removed from segments that begin at an in-silico cut before those segments are paired or written to `--segments-out`.

Importantly, stripping is performed on the **segment itself**, not specifically on R1 or R2.

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

The unique `--segments-out` file always represents the processed segment itself in native read orientation; the R2-specific reverse-complement option does not alter the unique segment output.

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

Pair identity therefore comes directly from matching FASTQ record names.

R1 and R2 do **not** use separate `/1` and `/2` suffixes.

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

The report includes information about both the original CiFi reads and the generated paired contacts.

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

When `--segments-out` is used, the command also reports the number of unique retained segments written.

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

Custom recognition sites can be supplied using:

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

`cifi contacts` accepts a name-sorted BAM containing mapped CiFi segments.

For aligned paired-contact filtering, use:

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
                                ▼
                          cifi digest
                                │
               ┌────────────────┴────────────────┐
               ▼                                 ▼
        R1/R2 paired FASTQ               unique segments
               │                         (--segments-out)
               ▼                                 │
            hifiasm                              ▼
      phasing / assembly                    minimap2
                                                  │
                                                  ▼
                                          samtools sort -n
                                                  │
                                                  ▼
                                           cifi contacts
                                                  │
                                                  ▼
                                                 BED
                                                  │
                                                  ▼
                                                YaHS
                                                  │
                                                  ▼
                                             scaffolds
```

The same CiFi digestion therefore provides:

```text
pair-expanded contacts → hifiasm
unique mapped segments → YaHS scaffolding
```

without requiring the same segment to be aligned repeatedly for every pairwise combination.

---

# Terminology

### CiFi read

The original long PacBio HiFi sequence produced from a proximity-ligated CiFi molecule.

### Segment

A sequence recovered from a CiFi read by in-silico digestion at the specified restriction-enzyme sites.

### Pair

Two segments from the same CiFi read emitted together as corresponding R1 and R2 FASTQ records.

### Contact

The pairwise relationship represented by two CiFi segments.

A CiFi read containing multiple segments therefore represents multiple pairwise contacts.

### Unique segment

A retained segment written once by:

```text
cifi digest --segments-out
```

and named:

```text
<read>__CIFI_SEG__<k>
```

Its mapped coordinate is used by `cifi contacts` to reconstruct the molecule's pairwise contacts.

---

# Citation

If you use `cifi-toolkit`, please cite:

**Single-library chromosome-scale diploid assemblies of vole genomes resolve a species-specific duplication implicated in pair bonding.**

*Cell Genomics.* 2026; 101336.

https://doi.org/10.1016/j.xgen.2026.101336

For the CiFi method, please also cite:

**CiFi: accurate long-read chromosome conformation capture with low-input requirements.**

*Nature Communications.* 2025.

https://doi.org/10.1038/s41467-025-66918-y

---

# License

MIT License

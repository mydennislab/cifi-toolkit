"""Command-line interface for CiFi - downstream processing of CiFi long reads."""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import click

from . import __version__


@click.group()
@click.version_option()
def main():
    """CiFi - toolkit for downstream processing of CiFi long reads.

    https://voles.dennislab.org
    """
    pass


def _validate_enzyme(ctx, param, value):
    """Validate enzyme name (only if provided)."""
    if value is None:
        return value
    from . import list_enzymes
    available = list_enzymes()
    if value not in available:
        raise click.BadParameter(
            f"Unknown enzyme '{value}'. Available: {', '.join(sorted(available))}"
        )
    return value


def _validate_site(ctx, param, value):
    """Validate custom enzyme site (IUPAC DNA sequence)."""
    if value is None:
        return value

    valid_bases = set("ACGTNRYWSMKBDHVacgtnrywsmkbdhv")
    for base in value:
        if base not in valid_bases:
            raise click.BadParameter(
                f"Invalid base '{base}' in site. Use IUPAC codes: A, C, G, T, N, R, Y, W, S, M, K, B, D, H, V"
            )
    return value.upper()


def _summarize(stats, fast_mode, histogram_bins=None, integer_valued=False):
    """Summary dict for a Statistics object, with an optional histogram.

    Histograms are binned in C++ so a report does not copy every recorded value
    into Python. For integer-valued data the bin count is capped at the number
    of distinct values, otherwise most bins come out structurally empty.
    """
    if stats.count() == 0:
        return None
    out = {
        "count": stats.count(),
        "min": stats.min(),
        "max": stats.max(),
        "mean": stats.mean(),
        "median": stats.median(),
    }
    if not fast_mode:
        out["q25"] = stats.percentile(0.25)
        out["q75"] = stats.percentile(0.75)
    if histogram_bins:
        n = histogram_bins
        if integer_valued:
            n = max(1, min(n, stats.max() - stats.min() + 1))
        edges, counts = stats.binned(n, integer_valued)
        if edges:
            out["histogram"] = {"bins": list(edges), "counts": list(counts)}
    return out


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "-e", "--enzyme", default=None, callback=_validate_enzyme,
    help="Restriction enzyme name. Use 'cifi enzymes' to list available enzymes."
)
@click.option(
    "--site", default=None, callback=_validate_site,
    help="Custom recognition site (supports IUPAC: N, R, Y, W, S, M, K, B, D, H, V)"
)
@click.option(
    "--cut-pos", "cut_offset", type=int, default=None,
    help="Cut position within site (0 to site length). Required with --site."
)
@click.option("-o", "--output-prefix", required=True, help="Output prefix for R1/R2 files")
@click.option(
    "-m", "--min-segments", default=3, show_default=True,
    help="Minimum segments required per read"
)
@click.option(
    "-l", "--min-segment-len", default=60, show_default=True,
    help="Minimum length of an emitted read (bp). Guaranteed for both R1 and R2."
)
@click.option(
    "--strip-overhang/--no-strip-overhang", default=True, show_default=True,
    help="Drop the 5' site remnant from segments that begin at a cut"
)
@click.option(
    "--revcomp-r2/--no-revcomp-r2", default=False, show_default=True,
    help="Reverse complement R2 (applied after stripping). Off by default: "
         "segments keep the orientation they were sequenced in."
)
@click.option(
    "--report/--no-report", default=True, show_default=True,
    help="Generate comprehensive HTML report"
)
@click.option(
    "--json/--no-json", "write_json", default=True, show_default=True,
    help="Write JSON statistics file"
)
@click.option(
    "--gzip", "gzip_output", is_flag=True, default=False,
    help="Compress output files with gzip"
)
@click.option(
    "--fast", "fast_mode", is_flag=True, default=False,
    help="Use streaming statistics (lower memory, approximate percentiles)"
)
def digest(input_file, enzyme, site, cut_offset, output_prefix, min_segments, min_segment_len,
           strip_overhang, revcomp_r2, report, write_json, gzip_output, fast_mode):
    """In-silico restriction digestion, generating paired-end FASTQ.

    Both mates of a pair share one read name, and -l bounds the length of the
    emitted reads, so R1 and R2 stay in step and neither drops below the cutoff.

    \b
    Examples:
        cifi digest reads.bam -e HindIII -o output
        cifi digest reads.fq.gz -e NlaIII -o output -m 5 --gzip
        cifi digest reads.bam --site GANTC --cut-pos 1 -o output
    """
    from . import get_enzyme_info, is_bam_file, process_reads, process_reads_custom

    # Validate enzyme/site options
    if enzyme is None and site is None:
        raise click.UsageError("Either --enzyme or --site is required.")

    if enzyme is not None and site is not None:
        raise click.UsageError("Use either --enzyme or --site, not both.")

    if site is not None and cut_offset is None:
        raise click.UsageError("--cut-pos is required when using --site.")

    if site is not None and cut_offset is not None:
        if cut_offset < 0 or cut_offset > len(site):
            raise click.UsageError(f"--cut-pos must be between 0 and {len(site)} (site length).")

    # Get enzyme info
    use_custom = site is not None
    if use_custom:
        enzyme_name = f"Custom({site})"
    else:
        site, cut_offset = get_enzyme_info(enzyme)
        enzyme_name = enzyme
    fmt = "BAM/SAM" if is_bam_file(input_file) else "FASTQ"

    # Determine if gzip output (explicit flag OR prefix ends with .gz)
    use_gzip = gzip_output or output_prefix.endswith('.gz')
    if output_prefix.endswith('.gz'):
        output_prefix = output_prefix[:-3]  # Remove .gz for prefix

    # Output paths
    ext = ".fastq.gz" if use_gzip else ".fastq"
    out_r1 = f"{output_prefix}_R1{ext}"
    out_r2 = f"{output_prefix}_R2{ext}"

    click.echo("CiFi Digestion Report")
    click.echo("-" * 60)
    click.echo(f"Input:       {input_file}")
    click.echo(f"Format:      {fmt}")
    click.echo(f"Enzyme:      {enzyme_name} ({site})")
    click.echo(f"Cut position: {cut_offset} (0-indexed)")
    click.echo(f"Min segments:   {min_segments}")
    click.echo(f"Min read len: {min_segment_len} bp (emitted R1/R2)")
    click.echo(f"Overhang:    {'stripped' if strip_overhang else 'kept'}")
    click.echo(f"R2 strand:   {'reverse complement' if revcomp_r2 else 'native'}")
    _trim = len(site) - cut_offset if 0 < cut_offset < len(site) else (len(site) if cut_offset == 0 else 0)
    if strip_overhang and _trim == 0:
        click.echo(f"Warning: --strip-overhang has no effect for {enzyme_name}: the cut "
                   f"leaves no 5' remnant on the downstream segment.", err=True)
    click.echo(f"Compression: {'gzip' if use_gzip else 'none'}")
    click.echo(f"Stats mode:  {'fast (approximate)' if fast_mode else 'exact'}")
    click.echo("-" * 60)
    click.echo("Processing...", nl=False)

    try:
        if use_custom:
            result = process_reads_custom(
                input_file, out_r1, out_r2, site, cut_offset, min_segments, min_segment_len,
                strip_overhang, use_gzip, fast_mode, revcomp_r2
            )
        else:
            result = process_reads(
                input_file, out_r1, out_r2, enzyme, min_segments, min_segment_len, strip_overhang,
                use_gzip, fast_mode, revcomp_r2
            )
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        sys.exit(1)

    click.echo(" Done!")
    click.echo("-" * 60)

    # Display results
    if result.total_bases_in > 0:
        click.echo("\nINPUT READS:")
        click.echo(f"  Reads:             {result.reads_in:>12,}")
        click.echo(f"  Total bases:       {result.total_bases_in:>12,}")
        click.echo(f"  GC content:        {100.0 * result.gc_bases_in / result.total_bases_in:>11.1f}%")
        click.echo(f"  Enzyme sites:      {result.total_sites:>12,}")
        click.echo(f"  Sites/read (all):  {result.total_sites / result.reads_in:>12.1f}")
        rl = result.read_length_stats
        if rl.count() > 0:
            click.echo(f"  Read length:       {rl.min():,} - {rl.max():,} bp "
                       f"(mean {rl.mean():,.0f}, median {rl.median():,.0f})")

    click.echo("\nRESULTS:")
    click.echo(f"  Reads processed:   {result.reads_in:>12,}")
    click.echo(f"  Reads passing:     {result.reads_out:>12,}")
    click.echo(f"  Reads filtered:    {result.reads_skipped:>12,}")

    if result.reads_in > 0:
        pass_rate = 100 * result.reads_out / result.reads_in
        click.echo(f"  Pass rate:         {pass_rate:>11.1f}%")

    click.echo(f"\n  Total segments:   {result.total_segments:>12,}")
    click.echo(f"  Total pairs:       {result.pairs_written:>12,}")

    if result.reads_out > 0:
        avg_segments = result.total_segments / result.reads_out
        avg_pairs = result.pairs_written / result.reads_out
        click.echo(f"  Avg segments/read:    {avg_segments:>12.1f}")
        click.echo(f"  Avg pairs/read:    {avg_pairs:>12.1f}")

    # Segment length statistics (using Statistics object)
    if result.segment_length_stats.count() > 0:
        segment_stats = result.segment_length_stats
        click.echo("\nSEGMENT LENGTHS:")
        click.echo(f"  Range:   {segment_stats.min():,} - {segment_stats.max():,} bp")
        click.echo(f"  Mean:    {segment_stats.mean():,.0f} bp")
        click.echo(f"  Median:  {segment_stats.median():,.0f} bp")
        if not fast_mode:
            click.echo(f"  IQR:     {segment_stats.percentile(0.25):,.0f} - {segment_stats.percentile(0.75):,.0f} bp")

    # Sites per read statistics (using Statistics object)
    if result.sites_per_read_stats.count() > 0:
        sites_stats = result.sites_per_read_stats
        click.echo(f"\nSITES PER READ ({enzyme_name}, passing reads):")
        click.echo(f"  Range:   {sites_stats.min()} - {sites_stats.max()}")
        click.echo(f"  Mean:    {sites_stats.mean():.1f}")
        click.echo(f"  Median:  {sites_stats.median():.0f}")

    bases_in_segments = result.segment_length_stats.sum()
    bases_out = result.bases_out_r1 + result.bases_out_r2
    if result.total_bases_in > 0:
        click.echo("\nYIELD:")
        click.echo(f"  Bases in:          {result.total_bases_in:>12,}")
        click.echo(f"  Bases in segments: {bases_in_segments:>12,} "
                   f"({100.0 * bases_in_segments / result.total_bases_in:.1f}% retained)")
        if bases_in_segments:
            click.echo(f"  Bases written:     {bases_out:>12,} "
                       f"({bases_out / bases_in_segments:.1f}x, segments reused across pairs)")
        click.echo(f"  Trimmed overhang:  {result.bases_trimmed_overhang:>12,}")
        click.echo(f"  Dropped (< {min_segment_len}bp):  {result.bases_dropped_short:>12,} "
                   f"in {result.segments_dropped_short:,} segments")
        click.echo(f"  In filtered reads: {result.bases_in_filtered_reads:>12,}")

    click.echo("-" * 60)
    click.echo("\nOUTPUT FILES:")
    click.echo(f"  {out_r1}")
    click.echo(f"  {out_r2}")

    # Build stats data for JSON and report
    stats_data = {
        "cifi_version": __version__,
        "timestamp": datetime.now().isoformat(),
        "input": {
            "file": os.path.basename(input_file),
            "path": os.path.abspath(input_file),
            "format": fmt,
        },
        "parameters": {
            "enzyme": enzyme_name,
            "enzyme_site": site,
            "cut_offset": cut_offset,
            "custom_enzyme": use_custom,
            "min_segments": min_segments,
            "min_segment_len": min_segment_len,
            "strip_overhang": strip_overhang,
            "revcomp_r2": revcomp_r2,
            "gzip_output": use_gzip,
            "fast_mode": fast_mode,
        },
        "results": {
            "reads_in": result.reads_in,
            "reads_out": result.reads_out,
            "reads_skipped": result.reads_skipped,
            "filtered_few_sites": result.filtered_few_sites,
            "filtered_short_segments": result.filtered_short_segments,
            "pairs_written": result.pairs_written,
            "total_segments": result.total_segments,
            "pass_rate": result.reads_out / result.reads_in if result.reads_in > 0 else 0,
            "avg_segments_per_read": result.total_segments / result.reads_out if result.reads_out > 0 else 0,
            "avg_pairs_per_read": result.pairs_written / result.reads_out if result.reads_out > 0 else 0,
        },
        "output": {
            "r1": out_r1,
            "r2": out_r2,
        },
    }

    # Segment length statistics from Statistics object
    if result.segment_length_stats.count() > 0:
        stats_data["segment_lengths"] = {
            "count": result.segment_length_stats.count(),
            "min": result.segment_length_stats.min(),
            "max": result.segment_length_stats.max(),
            "mean": result.segment_length_stats.mean(),
            "median": result.segment_length_stats.median(),
        }
        if not fast_mode:
            stats_data["segment_lengths"]["q25"] = result.segment_length_stats.percentile(0.25)
            stats_data["segment_lengths"]["q75"] = result.segment_length_stats.percentile(0.75)
            # Histogram from raw values (only available in exact mode)
    seg_hist = _summarize(result.segment_length_stats, fast_mode, 50)
    if seg_hist and "histogram" in seg_hist:
        stats_data["segment_length_histogram"] = seg_hist["histogram"]

    # Input profile: every read, including ones the filters later skipped
    bases_in = result.total_bases_in
    stats_data["input_reads"] = {
        "count": result.reads_in,
        "total_bases": bases_in,
        "gc_content": 100.0 * result.gc_bases_in / bases_in if bases_in else 0,
        "total_sites": result.total_sites,
        "mean_sites_per_read_all": result.total_sites / result.reads_in if result.reads_in else 0,
        "length": _summarize(result.read_length_stats, fast_mode, 50),
    }

    # Yield. Segments are emitted once per pair they take part in, so the bases
    # written exceed the unique segment bases; keep the two apart.
    bases_in_segments = result.segment_length_stats.sum()
    bases_out = result.bases_out_r1 + result.bases_out_r2
    stats_data["yield"] = {
        "bases_in": bases_in,
        "bases_in_segments": bases_in_segments,
        "bases_in_filtered_reads": result.bases_in_filtered_reads,
        "fraction_retained": bases_in_segments / bases_in if bases_in else 0,
        "bases_out_r1": result.bases_out_r1,
        "bases_out_r2": result.bases_out_r2,
        "bases_out_total": bases_out,
        "expansion_factor": bases_out / bases_in_segments if bases_in_segments else 0,
    }

    # Where the input that did not become segments went
    stats_data["filtering"] = {
        "reads_skipped_few_sites": result.filtered_few_sites,
        "reads_skipped_short_segments": result.filtered_short_segments,
        "segments_dropped_short": result.segments_dropped_short,
        "bases_dropped_short_segments": result.bases_dropped_short,
        "bases_trimmed_overhang": result.bases_trimmed_overhang,
        "bases_in_filtered_reads": result.bases_in_filtered_reads,
    }

    # Distributions behind the averages (passing reads only)
    for key, stats in (("segments_per_read", result.segments_per_read_stats),
                       ("pairs_per_read", result.pairs_per_read_stats)):
        summary = _summarize(stats, fast_mode, 50, integer_valued=True)
        if summary:
            stats_data[key] = summary

    # Sites per read statistics from Statistics object
    if result.sites_per_read_stats.count() > 0:
        stats_data["sites_per_read"] = {
            "count": result.sites_per_read_stats.count(),
            "min": result.sites_per_read_stats.min(),
            "max": result.sites_per_read_stats.max(),
            "mean": result.sites_per_read_stats.mean(),
            "median": result.sites_per_read_stats.median(),
        }
        sites_hist = _summarize(result.sites_per_read_stats, fast_mode, 50,
                                integer_valued=True)
        if sites_hist and "histogram" in sites_hist:
            stats_data["sites_per_read_histogram"] = sites_hist["histogram"]

    # Write JSON stats
    if write_json:
        json_file = f"{output_prefix}_stats.json"
        with open(json_file, "w") as f:
            json.dump(stats_data, f, indent=2)
        click.echo(f"  {json_file}")

    # Generate HTML report
    if report:
        try:
            from .report import generate_digest_report
            report_file = f"{output_prefix}_digestion_report.html"
            generate_digest_report(stats_data, report_file)
            click.echo(f"  {report_file}")
        except ImportError as e:
            click.echo(f"\nWarning: Report skipped (missing jinja2): {e}", err=True)
        except Exception as e:
            click.echo(f"\nWarning: Report failed: {e}", err=True)


@main.command("qc")
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "-e", "--enzyme", default=None, callback=_validate_enzyme,
    help="Restriction enzyme to analyze"
)
@click.option(
    "--site", default=None, callback=_validate_site,
    help="Custom recognition site (supports IUPAC: N, R, Y, W, S, M, K, B, D, H, V)"
)
@click.option(
    "--cut-pos", "cut_offset", type=int, default=None,
    help="Cut position within site (0 to site length). Required with --site."
)
@click.option(
    "-o", "--output", default=None,
    help="Output prefix for reports (default: input filename)"
)
@click.option(
    "-n", "--num-reads", default=10000, show_default=True,
    help="Number of reads to sample for QC (0 = all reads)"
)
@click.option(
    "--min-sites", default=2, show_default=True,
    help="Minimum sites for a read to be considered 'usable'"
)
@click.option(
    "--html/--no-html", default=True, show_default=True,
    help="Generate interactive HTML report"
)
@click.option(
    "--json/--no-json", "write_json", default=True, show_default=True,
    help="Generate JSON results file"
)
@click.option(
    "--pdf/--no-pdf", default=True, show_default=True,
    help="Generate PDF report (requires matplotlib)"
)
@click.option(
    "-q", "--quiet", is_flag=True,
    help="Suppress terminal output"
)
def qc(input_file, enzyme, site, cut_offset, output, num_reads, min_sites, html, write_json, pdf, quiet):
    """Sample reads and report enzyme site frequency, segment sizes, and estimated yield.

    \b
    Examples:
        cifi qc reads.bam -e HindIII -o qc_out
        cifi qc reads.bam -e NlaIII -n 50000 -o qc_out
        cifi qc reads.bam --site GANTC --cut-pos 1 -o qc_out
    """
    from . import get_enzyme_info, is_bam_file

    # Validate enzyme/site options
    if enzyme is None and site is None:
        raise click.UsageError("Either --enzyme or --site is required.")

    if enzyme is not None and site is not None:
        raise click.UsageError("Use either --enzyme or --site, not both.")

    if site is not None and cut_offset is None:
        raise click.UsageError("--cut-pos is required when using --site.")

    if site is not None and cut_offset is not None:
        if cut_offset < 0 or cut_offset > len(site):
            raise click.UsageError(f"--cut-pos must be between 0 and {len(site)} (site length).")

    # Determine output prefix
    if output is None:
        output = Path(input_file).stem

    # Get enzyme info
    use_custom = site is not None
    if use_custom:
        enzyme_name = f"Custom({site})"
    else:
        site, cut_offset = get_enzyme_info(enzyme)
        enzyme_name = enzyme

    fmt = "BAM/SAM" if is_bam_file(input_file) else "FASTQ"

    if not quiet:
        click.echo("CiFi QC Analysis")
        click.echo("-" * 60)
        click.echo(f"Input:       {input_file}")
        click.echo(f"Format:      {fmt}")
        click.echo(f"Enzyme:      {enzyme_name} ({site})")
        click.echo(f"Cut position: {cut_offset} (0-indexed)")
        click.echo(f"Sampling:    {'all reads' if num_reads == 0 else f'{num_reads:,} reads'}")
        click.echo("─" * 60)
        click.echo("Analyzing...", nl=False)

    # Run QC analysis
    try:
        qc_result = _run_qc_for_enzyme(input_file, site, cut_offset, num_reads, min_sites)
    except Exception as e:
        if not quiet:
            click.echo(f"\nError: {e}", err=True)
        sys.exit(1)

    if qc_result["reads_analyzed"] == 0:
        if not quiet:
            click.echo("\nError: No reads found in input file", err=True)
        sys.exit(1)

    if not quiet:
        click.echo(f" Done! ({qc_result['reads_analyzed']:,} reads)")
        click.echo("─" * 60)

        # Display summary
        click.echo("\nREAD STATISTICS:")
        click.echo(f"  Reads analyzed:    {qc_result['reads_analyzed']:>12,}")
        click.echo(f"  Total bases:       {qc_result['total_bases']:>12,} ({qc_result['total_bases']/1e9:.2f} Gb)")
        click.echo(f"  Avg read length:   {qc_result['avg_read_length']:>12,.0f} bp")
        click.echo(f"  Median length:     {qc_result['median_read_length']:>12,.0f} bp")
        click.echo(f"  GC content:        {qc_result['gc_content']:>11.1f}%")

        click.echo(f"\nENZYME ANALYSIS ({enzyme_name}):")
        click.echo(f"  Recognition site:  {site}")
        click.echo(f"  Total sites found: {qc_result['total_sites']:>12,}")
        click.echo(f"  Sites per read:    {qc_result['sites_per_read_mean']:>12.1f} (mean)")
        click.echo(f"                     {qc_result['sites_per_read_median']:>12.0f} (median)")

        click.echo("\nESTIMATED YIELD:")
        click.echo(f"  Reads with ≥{min_sites} sites: {qc_result['reads_passing']:>10,} ({qc_result['pass_rate']:.1f}%)")
        click.echo(f"  Est. segments:    {qc_result['est_total_segments']:>12,}")
        click.echo(f"  Est. pairs:        {qc_result['est_total_pairs']:>12,}")
        click.echo(f"  Avg segments/read:    {qc_result['avg_segments_per_read']:>12.1f}")
        click.echo(f"  Avg pairs/read:    {qc_result['avg_pairs_per_read']:>12.1f}")

        if qc_result.get('segment_size_mean'):
            click.echo("\nSEGMENT SIZE ESTIMATES:")
            click.echo(f"  Mean:              {qc_result['segment_size_mean']:>12,.0f} bp")
            click.echo(f"  Median:            {qc_result['segment_size_median']:>12,.0f} bp")

        click.echo("─" * 60)

    # Build full results for output
    results = {
        "cifi_version": __version__,
        "timestamp": datetime.now().isoformat(),
        "analysis_type": "qc",
        "input": {
            "file": os.path.basename(input_file),
            "path": os.path.abspath(input_file),
            "format": fmt,
        },
        "parameters": {
            "enzyme": enzyme_name,
            "enzyme_site": site,
            "cut_offset": cut_offset,
            "custom_enzyme": use_custom,
            "num_reads_sampled": num_reads,
            "min_sites_threshold": min_sites,
        },
        **qc_result,
    }

    # Create output directory for all QC files
    qc_dir = output
    Path(qc_dir).mkdir(parents=True, exist_ok=True)
    output_files = []

    # Write JSON
    if write_json:
        json_path = os.path.join(qc_dir, "qc.json")
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2)
        output_files.append(json_path)

    # Generate HTML report
    if html:
        try:
            from .report import generate_qc_report
            html_path = os.path.join(qc_dir, "qc.html")
            generate_qc_report(results, html_path)
            output_files.append(html_path)
        except ImportError as e:
            if not quiet:
                click.echo(f"Warning: HTML report skipped: {e}", err=True)
        except Exception as e:
            if not quiet:
                click.echo(f"Warning: HTML report failed: {e}", err=True)

    # Write TSV files
    try:
        from .report import write_qc_tsvs
        tsv_files = write_qc_tsvs(results, qc_dir)
        output_files.extend(tsv_files)
    except Exception as e:
        if not quiet:
            click.echo(f"Warning: TSV writing failed: {e}", err=True)

    # Generate plot PNGs
    try:
        from .report import generate_qc_plots
        plot_files = generate_qc_plots(results, qc_dir)
        output_files.extend(plot_files)
    except ImportError:
        pass  # matplotlib not available
    except Exception as e:
        if not quiet:
            click.echo(f"Warning: Plot generation failed: {e}", err=True)

    # Generate PDF report
    if pdf:
        try:
            from .report import generate_qc_pdf
            pdf_path = os.path.join(qc_dir, "qc.pdf")
            result_path = generate_qc_pdf(results, pdf_path)
            if result_path:
                output_files.append(result_path)
            elif not quiet:
                click.echo("Warning: PDF skipped (matplotlib not available)", err=True)
        except ImportError:
            if not quiet:
                click.echo("Warning: PDF skipped (matplotlib not available)", err=True)
        except Exception as e:
            if not quiet:
                click.echo(f"Warning: PDF generation failed: {e}", err=True)

    if not quiet and output_files:
        click.echo("\nOutput files:")
        for f in output_files:
            click.echo(f"  {f}")


def _run_qc_for_enzyme(input_file: str, site: str, cut_offset: int, num_reads: int, min_sites: int) -> dict:
    """Run QC analysis for a specific enzyme site using C++ backend.

    Args:
        input_file: Path to input FASTQ/BAM file
        site: Recognition sequence (supports IUPAC degenerate bases)
        cut_offset: Cut position within site
        num_reads: Number of reads to sample
        min_sites: Minimum sites for a read to be considered 'usable'

    Returns a dict with all QC metrics.
    """
    from . import run_qc_analysis_custom

    # Call C++ backend
    result = run_qc_analysis_custom(input_file, site, cut_offset, num_reads, min_sites)

    if result.reads_analyzed == 0:
        return {"reads_analyzed": 0}

    # Convert C++ result to dict format expected by the rest of the code
    return {
        "reads_analyzed": result.reads_analyzed,
        "total_bases": result.total_bases,
        "avg_read_length": result.avg_read_length,
        "median_read_length": result.median_read_length,
        "min_read_length": result.min_read_length,
        "max_read_length": result.max_read_length,
        "gc_content": result.gc_content,

        "total_sites": result.total_sites,
        "sites_per_read_mean": result.sites_per_read_mean,
        "sites_per_read_median": result.sites_per_read_median,
        "sites_per_read_min": result.sites_per_read_min,
        "sites_per_read_max": result.sites_per_read_max,

        "reads_passing": result.reads_passing,
        "pass_rate": result.pass_rate,
        "est_total_segments": result.est_total_segments,
        "est_total_pairs": result.est_total_pairs,
        "avg_segments_per_read": result.avg_segments_per_read,
        "avg_pairs_per_read": result.avg_pairs_per_read,

        "segment_size_mean": result.segment_size_mean,
        "segment_size_median": result.segment_size_median,
        "segment_size_min": result.segment_size_min,
        "segment_size_max": result.segment_size_max,

        # Histograms from C++ result
        "read_length_histogram": {
            "bins": list(result.read_length_hist_bins),
            "counts": list(result.read_length_hist_counts)
        } if result.read_length_hist_bins else {"bins": [], "counts": []},

        "sites_histogram": {
            "bins": list(result.sites_hist_bins),
            "counts": list(result.sites_hist_counts)
        } if result.sites_hist_bins else {"bins": [], "counts": []},

        "segment_size_histogram": {
            "bins": list(result.segment_size_hist_bins),
            "counts": list(result.segment_size_hist_counts)
        } if result.segment_size_hist_bins else None,
    }


def _make_histogram(data, num_bins):
    """Create histogram data for JSON output."""
    if not data:
        return {"bins": [], "counts": []}

    min_val, max_val = min(data), max(data)
    if min_val == max_val:
        return {"bins": [min_val], "counts": [len(data)]}

    bin_width = (max_val - min_val) / num_bins
    bins = [min_val + i * bin_width for i in range(num_bins + 1)]
    counts = [0] * num_bins

    for val in data:
        idx = min(int((val - min_val) / bin_width), num_bins - 1)
        counts[idx] += 1

    return {"bins": bins, "counts": counts}


@main.command()
def enzymes():
    """List available restriction enzymes."""
    from . import get_enzyme_info, list_enzymes

    click.echo("CiFi - Available Restriction Enzymes")
    click.echo("-" * 45)

    # Group by cutter type
    cutters = {"4": [], "6": [], "8": []}

    for name in list_enzymes():
        site, offset = get_enzyme_info(name)
        # Show cut position
        display = site[:offset] + "↓" + site[offset:]
        entry = (name, display, site, offset)

        if len(site) <= 4:
            cutters["4"].append(entry)
        elif len(site) <= 6:
            cutters["6"].append(entry)
        else:
            cutters["8"].append(entry)

    click.echo("\n4-CUTTERS (frequent cuts, many small segments):")
    click.echo(f"  {'Name':<10} {'Cut Site':<12} {'Notes'}")
    click.echo("  " + "─" * 40)
    for name, display, site, _ in sorted(cutters["4"]):
        notes = ""
        if name in ("DpnII", "MboI", "Sau3AI"):
            notes = "GATC cutters (common)"
        click.echo(f"  {name:<10} {display:<12} {notes}")

    click.echo("\n6-CUTTERS (moderate cuts, standard for CiFi):")
    click.echo(f"  {'Name':<10} {'Cut Site':<12} {'Notes'}")
    click.echo("  " + "─" * 40)
    for name, display, site, _ in sorted(cutters["6"]):
        notes = ""
        if name == "HindIII":
            notes = "common CiFi enzyme"
        click.echo(f"  {name:<10} {display:<12} {notes}")

    if cutters["8"]:
        click.echo("\n8-CUTTERS (rare cuts, few segments):")
        click.echo(f"  {'Name':<10} {'Cut Site':<12}")
        click.echo("  " + "─" * 40)
        for name, display, _, _ in sorted(cutters["8"]):
            click.echo(f"  {name:<10} {display:<12}")

    click.echo("\n↓ indicates the cut position in the recognition sequence")
    click.echo("\nUsage: cifi qc reads.bam -e <ENZYME>")


@main.command("filter")
@click.argument("input_bam", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output BAM file path")
@click.option("-q", "--mapq", default=30, show_default=True,
              help="Minimum MAPQ threshold")
@click.option("-t", "--threads", default=4, show_default=True,
              help="Number of threads for BAM I/O")
@click.option("--report/--no-report", default=True, show_default=True,
              help="Generate HTML report")
@click.option("--json/--no-json", "write_json", default=True, show_default=True,
              help="Write JSON statistics file")
@click.option("--quiet", is_flag=True, help="Suppress terminal output")
def filter_cmd(input_bam, output, mapq, threads, report, write_json, quiet):
    """Filter aligned paired-end BAM by mapping quality.

    \b
    Example:
        cifi filter aligned.bam -o filtered.bam -q 30
    """
    from . import filter_bam

    result = filter_bam(input_bam, output, mapq, threads)

    pass_rate = result.passed_pairs / result.total_pairs * 100 if result.total_pairs > 0 else 0

    if not quiet:
        click.echo("\nFilter Summary")
        click.echo(f"{'─' * 40}")
        click.echo(f"Total reads:      {result.total_reads:,}")
        click.echo(f"Total pairs:      {result.total_pairs:,}")
        click.echo(f"Passed pairs:     {result.passed_pairs:,} ({pass_rate:.1f}%)")
        click.echo(f"{'─' * 40}")
        click.echo(f"Failed - unpaired:  {result.failed_unpaired:,}")
        click.echo(f"Failed - unmapped:  {result.failed_unmapped:,}")
        click.echo(f"Failed - low MAPQ:  {result.failed_mapq:,}")
        if result.failed_mate_not_found > 0:
            click.echo(f"Failed - mate not found: {result.failed_mate_not_found:,}")

    # Build stats data
    output_prefix = output.replace('.bam', '').replace('.BAM', '')
    stats_data = {
        "cifi_version": __version__,
        "timestamp": datetime.now().isoformat(),
        "input": {
            "file": os.path.basename(input_bam),
            "path": os.path.abspath(input_bam),
        },
        "parameters": {
            "mapq_threshold": mapq,
            "threads": threads,
        },
        "results": {
            "total_reads": result.total_reads,
            "total_pairs": result.total_pairs,
            "passed_pairs": result.passed_pairs,
            "pass_rate": pass_rate,
            "failed_unpaired": result.failed_unpaired,
            "failed_unmapped": result.failed_unmapped,
            "failed_mapq": result.failed_mapq,
            "failed_mate_not_found": result.failed_mate_not_found,
        },
        "output": {
            "file": output,
        },
    }

    # Write JSON
    output_files = [output]
    if write_json:
        json_file = f"{output_prefix}_filter_stats.json"
        with open(json_file, "w") as f:
            json.dump(stats_data, f, indent=2)
        output_files.append(json_file)

    # Generate HTML report
    if report:
        try:
            from .report import generate_filter_report
            report_file = f"{output_prefix}_filter_report.html"
            generate_filter_report(stats_data, report_file)
            output_files.append(report_file)
        except ImportError as e:
            if not quiet:
                click.echo(f"\nWarning: Report skipped (missing jinja2): {e}", err=True)
        except Exception as e:
            if not quiet:
                click.echo(f"\nWarning: Report failed: {e}", err=True)

    if not quiet:
        click.echo("\nOutput files:")
        for f in output_files:
            click.echo(f"  {f}")


if __name__ == "__main__":
    main()

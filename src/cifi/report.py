"""CiFi Report Generator - Generates HTML, PDF, TSV reports and plot PNGs."""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    from jinja2 import Environment, FileSystemLoader
    HAS_JINJA2 = True
except ImportError:
    HAS_JINJA2 = False


def _format_number(value: Any) -> str:
    """Format number with commas."""
    try:
        return f"{int(value):,}" if isinstance(value, (int, float)) else str(value)
    except (ValueError, TypeError):
        return str(value)


def _histogram_to_plot_data(hist_data: Optional[Dict]) -> Optional[Dict]:
    """Convert histogram data to format for Plotly."""
    if not hist_data or not hist_data.get("bins"):
        return None
    bins, counts = hist_data["bins"], hist_data["counts"]
    if len(bins) > 1:
        bin_centers = [b + (bins[1] - bins[0]) / 2 for b in bins[:-1]]
    else:
        bin_centers = bins
    return {"bins": bin_centers, "counts": counts}


def _format_enzyme_site(site: str, cut_offset: int) -> str:
    """Format enzyme recognition site with cut position marker."""
    if 0 <= cut_offset <= len(site):
        return site[:cut_offset] + "|" + site[cut_offset:]
    return site


def _load_template() -> Any:
    """Load the HTML template from file."""
    if not HAS_JINJA2:
        raise ImportError("Jinja2 required for report generation: pip install jinja2")

    template_dir = Path(__file__).parent / "templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    return env.get_template("report.html")


def generate_qc_report(results: Dict[str, Any], output_path: str) -> str:
    """Generate HTML QC report.
    """
    template = _load_template()

    enzyme = results["parameters"]["enzyme"]
    site = results.get("enzyme_site", results["parameters"].get("enzyme_site", ""))
    cut_offset = results.get("cut_offset", results["parameters"].get("cut_offset", 0))

    # Build template data
    template_data = {
        "report_type": "QC",
        "enzyme": enzyme,
        "input_file": results["input"]["file"],
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),

        # Summary metrics
        "summary_metrics": [
            {"value": _format_number(results["reads_analyzed"]), "label": "Reads Analyzed"},
            {"value": f'{results["total_bases"] / 1e9:.2f} Gb', "label": "Total Bases"},
            {"value": _format_number(int(results["avg_read_length"])), "label": "Avg Read Length (bp)"},
            {"value": f'{results["gc_content"]:.1f}%', "label": "GC Content"},
        ],

        # Enzyme analysis table
        "enzyme_table": [
            ["Enzyme", enzyme],
            ["Recognition Site", _format_enzyme_site(site, cut_offset) if site else "N/A"],
            ["Total Sites Found", _format_number(results.get("total_sites", 0))],
            ["Sites per Read (mean)", f'{results.get("sites_per_read_mean", 0):.1f}'],
            ["Sites per Read (median)", f'{results.get("sites_per_read_median", 0):.0f}'],
            ["Pass Rate", f'{results.get("pass_rate", 0):.1f}%'],
            ["Estimated Total Pairs", _format_number(results.get("est_total_pairs", 0))],
        ],

        # Histograms
        "read_length_histogram": _histogram_to_plot_data(results.get("read_length_histogram")),
        "sites_histogram": _histogram_to_plot_data(results.get("sites_histogram")),
        "segment_size_histogram": _histogram_to_plot_data(results.get("segment_size_histogram")),
    }

    html = template.render(**template_data)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    return str(output_path)


def _param(params, new_key, old_key, default="n/a"):
    """Read a parameter, tolerating stats files written before the rename."""
    if new_key in params:
        return params[new_key]
    return params.get(old_key, default)


def generate_digest_report(stats_data: Dict[str, Any], output_path: str) -> str:
    """Generate HTML report for digest command results.
    """
    template = _load_template()

    results = stats_data["results"]
    params = stats_data["parameters"]

    # Calculate pass rate
    pass_rate = results["pass_rate"] * 100 if results["pass_rate"] <= 1 else results["pass_rate"]

    # Build filter table if applicable
    reads_in = results["reads_in"]
    filtered_few_sites = results.get("filtered_few_sites", 0)
    filtered_short_segments = results.get(
        "filtered_short_segments", results.get("filtered_short_frags", 0))

    filter_table = []
    if reads_in > 0 and (filtered_few_sites > 0 or filtered_short_segments > 0):
        if filtered_few_sites > 0:
            pct = 100 * filtered_few_sites / reads_in
            filter_table.append([
                f'Too few sites (< {_param(params, "min_segments", "min_fragments")} segments)',
                _format_number(filtered_few_sites),
                f'{pct:.1f}%'
            ])
        if filtered_short_segments > 0:
            pct = 100 * filtered_short_segments / reads_in
            filter_table.append([
                f'Segments too short (< {_param(params, "min_segment_len", "min_frag_len")} bp)',
                _format_number(filtered_short_segments),
                f'{pct:.1f}%'
            ])

    # Build segment stats table if available
    segment_stats = stats_data.get("segment_lengths") or stats_data.get("fragment_lengths")
    segment_stats_table = None
    if segment_stats:
        segment_stats_table = [
            ["Total Segments", _format_number(
                results.get("total_segments", results.get("total_fragments", 0)))],
            ["Mean Length", f'{segment_stats["mean"]:.0f} bp'],
            ["Median Length", f'{segment_stats["median"]:.0f} bp'],
            ["Range", f'{segment_stats["min"]:,} - {segment_stats["max"]:,} bp'],
        ]

    # Input profile: describes what was fed in, including skipped reads
    inp = stats_data.get("input_reads") or {}
    input_table = None
    if inp.get("total_bases"):
        length = inp.get("length") or {}
        input_table = [
            ["Reads", _format_number(inp["count"])],
            ["Total Bases", _format_number(inp["total_bases"])],
            ["GC Content", f'{inp["gc_content"]:.1f}%'],
            ["Enzyme Sites", _format_number(inp["total_sites"])],
            ["Mean Sites / Read (all)",
             f'{_param(inp, "mean_sites_per_read_all", "mean_sites_per_read", 0):.1f}'],
        ]
        if length:
            input_table += [
                ["Read Length (mean)", f'{length["mean"]:,.0f} bp'],
                ["Read Length (median)", f'{length["median"]:,.0f} bp'],
                ["Read Length (range)", f'{length["min"]:,} - {length["max"]:,} bp'],
            ]

    # Yield: unique segment bases vs bases written, which differ because each
    # segment is emitted once per pair it takes part in
    y = stats_data.get("yield") or {}
    yield_table = None
    if y.get("bases_in"):
        yield_table = [
            ["Bases In", _format_number(y["bases_in"])],
            ["Bases In Segments", f'{_format_number(y["bases_in_segments"])} '
                                  f'({100 * y["fraction_retained"]:.1f}% retained)'],
            ["Bases Written (R1)", _format_number(y["bases_out_r1"])],
            ["Bases Written (R2)", _format_number(y["bases_out_r2"])],
            ["Expansion Factor", f'{y["expansion_factor"]:.1f}x'],
        ]
        if "bases_in_filtered_reads" in y:
            yield_table.insert(2, ["In Filtered Reads",
                                   _format_number(y["bases_in_filtered_reads"])])
    f = stats_data.get("filtering") or {}
    if f:
        yield_table = (yield_table or []) + [
            ["Bases Trimmed (overhang)", _format_number(f["bases_trimmed_overhang"])],
            ["Bases Dropped (short segments)",
             f'{_format_number(f["bases_dropped_short_segments"])} '
             f'in {_format_number(f["segments_dropped_short"])} segments'],
        ]

    # Build template data
    template_data = {
        "report_type": "Digestion",
        "enzyme": params["enzyme"],
        "input_file": stats_data["input"]["file"],
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),

        # Parameter table
        "param_table": [
            ["Enzyme", params["enzyme"]],
            ["Recognition Site", params.get("enzyme_site", "N/A")],
            ["Min Segments", str(_param(params, "min_segments", "min_fragments"))],
            ["Min Segment Length", f'{_param(params, "min_segment_len", "min_frag_len")} bp'],
        ],

        # Summary metrics
        "summary_metrics": [
            {"value": _format_number(results["reads_in"]), "label": "Reads In"},
            {"value": _format_number(results["reads_out"]), "label": "Reads Passing"},
            {"value": _format_number(results["pairs_written"]), "label": "Pairs Written"},
            {"value": f'{pass_rate:.1f}%', "label": "Pass Rate"},
        ],

        # Filter table (may be empty)
        "filter_table": filter_table if filter_table else None,
        "filter_caption": f'Total filtered: {_format_number(results["reads_skipped"])} reads' if filter_table else "",

        # Segment stats table (may be None)
        "segment_stats_table": segment_stats_table,

        # Input profile and yield (may be None on older stats files)
        "input_table": input_table,
        "yield_table": yield_table,

        # Histograms
        "segment_length_histogram": _histogram_to_plot_data(stats_data.get("segment_length_histogram")
                                                or stats_data.get("fragment_length_histogram")),
        "sites_per_read_histogram": _histogram_to_plot_data(stats_data.get("sites_per_read_histogram")),
        "read_length_histogram": _histogram_to_plot_data((inp.get("length") or {}).get("histogram")),
        "segments_per_read_histogram": _histogram_to_plot_data(
            (stats_data.get("segments_per_read") or {}).get("histogram")),
        "pairs_per_read_histogram": _histogram_to_plot_data(
            (stats_data.get("pairs_per_read") or {}).get("histogram")),

        # Output files
        "output_files": [stats_data["output"]["r1"], stats_data["output"]["r2"]],
    }

    html = template.render(**template_data)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    return str(output_path)


def generate_filter_report(stats_data: Dict[str, Any], output_path: str) -> str:
    """Generate HTML report for filter command results.

    Args:
        stats_data: Dictionary containing filter statistics
        output_path: Path to write the HTML report

    Returns:
        The output path of the generated report
    """
    template = _load_template()

    results = stats_data["results"]
    params = stats_data["parameters"]

    # Build filter breakdown table
    filter_table = []
    total_failed = results["total_pairs"] - results["passed_pairs"]
    if total_failed > 0:
        if results["failed_unpaired"] > 0:
            pct = 100 * results["failed_unpaired"] / results["total_pairs"]
            filter_table.append(["Unpaired reads", _format_number(results["failed_unpaired"]), f"{pct:.1f}%"])
        if results["failed_unmapped"] > 0:
            pct = 100 * results["failed_unmapped"] / results["total_pairs"]
            filter_table.append(["Unmapped", _format_number(results["failed_unmapped"]), f"{pct:.1f}%"])
        if results["failed_mapq"] > 0:
            pct = 100 * results["failed_mapq"] / results["total_pairs"]
            filter_table.append([f"MAPQ < {params['mapq_threshold']}", _format_number(results["failed_mapq"]), f"{pct:.1f}%"])
        if results.get("failed_mate_not_found", 0) > 0:
            pct = 100 * results["failed_mate_not_found"] / results["total_pairs"]
            filter_table.append(["Mate not found", _format_number(results["failed_mate_not_found"]), f"{pct:.1f}%"])

    template_data = {
        "report_type": "Filter",
        "enzyme": f"MAPQ ≥ {params['mapq_threshold']}",
        "input_file": stats_data["input"]["file"],
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),

        # Parameter table
        "param_table": [
            ["MAPQ Threshold", str(params["mapq_threshold"])],
            ["Threads", str(params["threads"])],
        ],

        # Summary metrics
        "summary_metrics": [
            {"value": _format_number(results["total_reads"]), "label": "Total Reads"},
            {"value": _format_number(results["total_pairs"]), "label": "Total Pairs"},
            {"value": _format_number(results["passed_pairs"]), "label": "Passed Pairs"},
            {"value": f'{results["pass_rate"]:.1f}%', "label": "Pass Rate"},
        ],

        # Filter breakdown table
        "filter_table": filter_table if filter_table else None,
        "filter_caption": f'Total filtered: {_format_number(total_failed)} pairs' if filter_table else "",

        # No histograms for filter
        "segment_length_histogram": None,
        "sites_per_read_histogram": None,

        # Output files
        "output_files": [stats_data["output"]["file"]],
    }

    html = template.render(**template_data)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    return str(output_path)


def write_qc_tsvs(results: Dict[str, Any], output_dir: str) -> List[str]:
    """Write QC results as TSV files to a subdirectory.

    Args:
        results: Dictionary containing QC analysis results
        output_dir: Directory to write TSV files into

    Returns:
        List of written file paths
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = []

    # summary.tsv — key-value pairs for scalar metrics
    summary_path = out / "summary.tsv"
    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        for key in [
            "reads_analyzed", "total_bases", "avg_read_length", "median_read_length",
            "min_read_length", "max_read_length", "gc_content",
            "total_sites", "sites_per_read_mean", "sites_per_read_median",
            "reads_passing", "pass_rate", "est_total_segments", "est_total_pairs",
            "avg_segments_per_read", "avg_pairs_per_read",
            "segment_size_mean", "segment_size_median", "segment_size_min", "segment_size_max",
        ]:
            val = results.get(key)
            if val is not None:
                f.write(f"{key}\t{val}\n")
    written.append(str(summary_path))

    # enzyme_analysis.tsv
    params = results.get("parameters", {})
    enzyme_path = out / "enzyme_analysis.tsv"
    with open(enzyme_path, "w") as f:
        f.write("field\tvalue\n")
        f.write(f"enzyme\t{params.get('enzyme', '')}\n")
        f.write(f"recognition_site\t{params.get('enzyme_site', '')}\n")
        f.write(f"cut_offset\t{params.get('cut_offset', '')}\n")
        f.write(f"total_sites\t{results.get('total_sites', 0)}\n")
        f.write(f"sites_per_read_mean\t{results.get('sites_per_read_mean', 0)}\n")
        f.write(f"sites_per_read_median\t{results.get('sites_per_read_median', 0)}\n")
        f.write(f"pass_rate\t{results.get('pass_rate', 0)}\n")
        f.write(f"est_total_pairs\t{results.get('est_total_pairs', 0)}\n")
    written.append(str(enzyme_path))

    # Histogram TSVs
    hist_specs = [
        ("read_length_histogram", "read_length_histogram.tsv"),
        ("sites_histogram", "sites_per_read_histogram.tsv"),
        ("segment_size_histogram", "segment_size_histogram.tsv"),
    ]
    for hist_key, filename in hist_specs:
        hist = results.get(hist_key)
        if not hist or not hist.get("bins"):
            continue
        bins = hist["bins"]
        counts = hist["counts"]
        path = out / filename
        with open(path, "w") as f:
            f.write("bin_start\tbin_end\tcount\n")
            for i in range(len(counts)):
                bin_start = bins[i]
                bin_end = bins[i + 1] if i + 1 < len(bins) else bins[i]
                f.write(f"{bin_start}\t{bin_end}\t{counts[i]}\n")
        written.append(str(path))

    return written


def generate_qc_plots(results: Dict[str, Any], output_dir: str) -> List[str]:
    """Generate 300 DPI plot PNGs for QC results.

    Args:
        results: Dictionary containing QC analysis results
        output_dir: Directory to write PNG files into

    Returns:
        List of written file paths, empty if matplotlib unavailable
    """
    if not HAS_MATPLOTLIB:
        return []

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = []

    plot_specs = [
        ("read_length_histogram", "read_length_distribution.png",
         "Read Length Distribution", "Read Length (bp)", "Count", "#3498db"),
        ("sites_histogram", "sites_per_read_distribution.png",
         "Sites per Read Distribution", "Sites per Read", "Count", "#27ae60"),
        ("segment_size_histogram", "segment_size_distribution.png",
         "Segment Size Distribution", "Segment Size (bp)", "Count", "#9b59b6"),
    ]

    for hist_key, filename, title, xlabel, ylabel, color in plot_specs:
        hist = results.get(hist_key)
        if not hist or not hist.get("bins") or not hist.get("counts"):
            continue

        bins = hist["bins"]
        counts = hist["counts"]

        if len(bins) > 1:
            width = bins[1] - bins[0]
            centers = [b + width / 2 for b in bins[:-1]]
        else:
            width = 1
            centers = bins

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(centers, counts, width=width * 0.9, color=color, edgecolor="white", linewidth=0.5)
        ax.set_title(title, fontsize=14, fontfamily="sans-serif", pad=12)
        ax.set_xlabel(xlabel, fontsize=12, fontfamily="sans-serif")
        ax.set_ylabel(ylabel, fontsize=12, fontfamily="sans-serif")
        ax.tick_params(labelsize=10)
        ax.set_facecolor("white")
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.patch.set_facecolor("white")

        path = out / filename
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        written.append(str(path))

    return written


def generate_qc_pdf(results: Dict[str, Any], output_path: str) -> Optional[str]:
    """Generate a multi-page PDF report for QC results.
    """
    if not HAS_MATPLOTLIB:
        return None

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    params = results.get("parameters", {})
    enzyme = params.get("enzyme", "")
    site = params.get("enzyme_site", "")
    cut_offset = params.get("cut_offset", 0)

    with PdfPages(output_path) as pdf:
        # Page 1: Title + Summary + Enzyme table
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis("off")

        # Title
        ax.text(0.5, 0.95, "CiFi QC Report", ha="center", va="top",
                fontsize=20, fontweight="bold", fontfamily="sans-serif")
        ax.text(0.5, 0.91, f"Enzyme: {enzyme}", ha="center", va="top",
                fontsize=14, fontfamily="sans-serif", color="#555")
        ax.text(0.5, 0.88, f"Input: {results.get('input', {}).get('file', '')}",
                ha="center", va="top", fontsize=11, fontfamily="sans-serif", color="#777")
        ax.text(0.5, 0.855, datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                ha="center", va="top", fontsize=10, fontfamily="sans-serif", color="#999")

        # Summary metrics table
        summary_data = [
            ["Reads Analyzed", _format_number(results.get("reads_analyzed", 0))],
            ["Total Bases", f'{results.get("total_bases", 0) / 1e9:.2f} Gb'],
            ["Avg Read Length", f'{results.get("avg_read_length", 0):,.0f} bp'],
            ["Median Read Length", f'{results.get("median_read_length", 0):,.0f} bp'],
            ["GC Content", f'{results.get("gc_content", 0):.1f}%'],
        ]

        tbl1 = ax.table(
            cellText=summary_data,
            colLabels=["Metric", "Value"],
            loc="center",
            bbox=[0.15, 0.55, 0.7, 0.25],
        )
        tbl1.auto_set_font_size(False)
        tbl1.set_fontsize(11)
        for key, cell in tbl1.get_celld().items():
            cell.set_edgecolor("#ddd")
            if key[0] == 0:
                cell.set_facecolor("#f0f0f0")
                cell.set_text_props(fontweight="bold")

        # Enzyme analysis table
        enzyme_data = [
            ["Enzyme", enzyme],
            ["Recognition Site", _format_enzyme_site(site, cut_offset) if site else "N/A"],
            ["Total Sites Found", _format_number(results.get("total_sites", 0))],
            ["Sites per Read (mean)", f'{results.get("sites_per_read_mean", 0):.1f}'],
            ["Sites per Read (median)", f'{results.get("sites_per_read_median", 0):.0f}'],
            ["Pass Rate", f'{results.get("pass_rate", 0):.1f}%'],
            ["Est. Total Pairs", _format_number(results.get("est_total_pairs", 0))],
        ]

        tbl2 = ax.table(
            cellText=enzyme_data,
            colLabels=["Enzyme Analysis", "Value"],
            loc="center",
            bbox=[0.15, 0.15, 0.7, 0.30],
        )
        tbl2.auto_set_font_size(False)
        tbl2.set_fontsize(11)
        for key, cell in tbl2.get_celld().items():
            cell.set_edgecolor("#ddd")
            if key[0] == 0:
                cell.set_facecolor("#f0f0f0")
                cell.set_text_props(fontweight="bold")

        pdf.savefig(fig, dpi=300)
        plt.close(fig)

        # Plot pages
        plot_specs = [
            ("read_length_histogram",
             "Read Length Distribution", "Read Length (bp)", "Count", "#3498db"),
            ("sites_histogram",
             "Sites per Read Distribution", "Sites per Read", "Count", "#27ae60"),
            ("segment_size_histogram",
             "Segment Size Distribution", "Segment Size (bp)", "Count", "#9b59b6"),
        ]

        for hist_key, title, xlabel, ylabel, color in plot_specs:
            hist = results.get(hist_key)
            if not hist or not hist.get("bins") or not hist.get("counts"):
                continue

            bins = hist["bins"]
            counts = hist["counts"]

            if len(bins) > 1:
                width = bins[1] - bins[0]
                centers = [b + width / 2 for b in bins[:-1]]
            else:
                width = 1
                centers = bins

            fig, ax = plt.subplots(figsize=(10, 7))
            ax.bar(centers, counts, width=width * 0.9, color=color,
                   edgecolor="white", linewidth=0.5)
            ax.set_title(title, fontsize=16, fontweight="bold",
                         fontfamily="sans-serif", pad=15)
            ax.set_xlabel(xlabel, fontsize=13, fontfamily="sans-serif")
            ax.set_ylabel(ylabel, fontsize=13, fontfamily="sans-serif")
            ax.tick_params(labelsize=11)
            ax.set_facecolor("white")
            ax.grid(axis="y", alpha=0.3, linewidth=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            fig.patch.set_facecolor("white")

            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)

    return str(output_path)

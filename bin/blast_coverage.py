#!/usr/bin/env python3
"""
BLAST Coverage Analysis Tool
Analyzes BLAST tabular output to calculate reference coverage and identity statistics.
Outputs MultiQC-compatible results for multiple samples.
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from textwrap import dedent

import typer
import numpy as np
import pandas as pd
from rich.logging import RichHandler


VERSION = "1.0.0"

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

app = typer.Typer()


OUTPUT_COLUMNS = [
    ("sample", "Sample"),
    ("reference", "Reference"),
    ("num_alignments", "# Alignments"),
    ("reference_length", "Reference Length"),
    ("total_alignment_length", "Total Alignment Length"),
    ("reference_coverage", "Reference Coverage"),
    ("avg_identity_covered", "Avg Identity Covered"),
    ("total_aligned_length", "Total Aligned Length"),
    ("reference_title", "Reference Title"),
]

def init_logging(verbose: bool) -> None:
    from rich.traceback import install
    install(show_locals=True, width=120, word_wrap=True)

    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG if verbose else logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


class BLASTCoverageAnalyzer:
    """Main class for BLAST coverage analysis."""

    def __init__(
        self,
        blast_file: str,
        output_dir: str = ".",
        top_n: int = 1,
    ):
        self.blast_file = blast_file
        self.output_dir = Path(output_dir)
        self.top_n = top_n
        self.subject_sequences = {}

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def parse_blast_output(self) -> pd.DataFrame:
        """Parse BLAST tabular output into a pandas DataFrame."""
        columns = [
            "qaccver",
            "saccver",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
            "stitle",
        ]

        try:
            df = pd.read_csv(self.blast_file, sep="\t", names=columns, header=None)
            logger.info(f"Loaded {len(df)} BLAST alignments from {self.blast_file}")
            return df
        except Exception as e:
            logger.error(f"Error parsing BLAST file {self.blast_file}: {e}")
            raise

    def calculate_reference_coverage(self, df: pd.DataFrame) -> Dict:
        """Calculate coverage and identity statistics for references."""
        results = {}

        # Group by query
        for sample, sample_group in df.groupby("qaccver"):
            # Get top N references by total bitscore
            ref_scores = (
                sample_group.groupby("saccver")["bitscore"]
                .sum()
                .sort_values(ascending=False)
            )
            top_refs = ref_scores[:self.top_n].index

            results[sample] = {}

            for ref in top_refs:
                # Get all alignments for this query-reference pair
                alignments = sample_group[sample_group["saccver"] == ref].copy()

                # Calculate coverage and identity statistics
                stats = self._calculate_alignment_stats(alignments, ref, sample)
                results[sample][ref] = stats

        return results

    def _calculate_alignment_stats(
        self, alignments: pd.DataFrame, ref_id: str, sample: str
    ) -> Dict:
        """Calculate detailed alignment statistics."""
        if len(alignments) == 0:
            return self._empty_stats(ref_id)

        # Get reference length from alignment 'slen' column
        ref_length = self._get_reference_length(ref_id, alignments)
        ref_title = alignments[alignments["saccver"] == ref_id]["stitle"].values[0]

        # Create coverage map and calculate statistics
        coverage_map = np.zeros(ref_length, dtype=bool)
        identity_map = np.zeros(ref_length, dtype=float)
        identity_counts = np.zeros(ref_length, dtype=int)

        regions = []
        weighted_identity_sum = 0

        for _, alignment in alignments.iterrows():
            start = (
                min(alignment["sstart"], alignment["send"]) - 1
            )  # Convert to 0-based
            end = max(alignment["sstart"], alignment["send"])
            length = alignment["length"]
            identity = alignment["pident"]
            identity = identity / 100.0

            regions.append((start + 1, end))  # Keep 1-based for output

            # Mark coverage
            coverage_map[start:end] = True

            # Track identity (weighted by alignment length)
            identity_map[start:end] += identity * length
            identity_counts[start:end] += length

            weighted_identity_sum += identity * length

        # Calculate final statistics
        covered_bases = np.sum(coverage_map)
        coverage_proportion = covered_bases / ref_length if ref_length > 0 else 0.0

        # Calculate average identity for covered regions only
        covered_positions = identity_counts > 0
        if np.sum(covered_positions) > 0:
            avg_identity_covered = np.sum(identity_map[covered_positions]) / np.sum(
                identity_counts[covered_positions]
            )
        else:
            avg_identity_covered = 0.0


        return {
            "reference_id": ref_id,
            "sample": sample,
            "total_alignments": len(alignments),
            "reference_length": ref_length,
            "total_alignment_length": int(covered_bases),
            "coverage_proportion": float(coverage_proportion),
            "avg_identity_covered_regions": float(avg_identity_covered),
            "reference_title": ref_title,
        }

    def _get_reference_length(self, ref_id: str, alignments: pd.DataFrame) -> int:
        """Get reference length from alignment 'slen' column."""
        return alignments[alignments["saccver"] == ref_id]["slen"].values[0]

    def _merge_overlapping_regions(
        self, regions: List[Tuple[int, int]]
    ) -> List[Tuple[int, int]]:
        """Merge overlapping genomic regions."""
        if not regions:
            return []

        sorted_regions = sorted(regions)
        merged = [sorted_regions[0]]

        for current_start, current_end in sorted_regions[1:]:
            last_start, last_end = merged[-1]

            if current_start <= last_end + 1:
                merged[-1] = (last_start, max(last_end, current_end))
            else:
                merged.append((current_start, current_end))

        return merged

    def _empty_stats(self, ref_id: str) -> Dict:
        """Return empty statistics dictionary."""
        return {
            "reference_id": ref_id,
            "sample": "",
            "total_alignments": 0,
            "reference_length": 0,
            "total_alignment_length": 0,
            "coverage_proportion": 0.0,
            "avg_identity_covered_regions": 0.0,
            "reference_title": "",
        }

    def create_summary_dataframe(self, coverage_results: Dict) -> pd.DataFrame:
        """Convert coverage results to a summary DataFrame."""
        rows = []
        for sample, references in coverage_results.items():
            for ref_id, stats in references.items():
                row = {
                    "sample": sample,
                    "reference": stats["reference_id"],
                    "num_alignments": stats["total_alignments"],
                    "reference_length": stats["reference_length"],
                    "total_alignment_length": stats["total_alignment_length"],
                    "reference_coverage": stats["coverage_proportion"],
                    "avg_identity_covered": stats["avg_identity_covered_regions"],
                    "reference_title": stats["reference_title"],
                }

                rows.append(row)

        return pd.DataFrame(rows)


class MultiSampleBLASTAnalyzer:
    """Analyze multiple BLAST files from a directory."""
    
    def __init__(self, input_dir: Path, output_dir: Path, top_n: int = 1):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.top_n = top_n
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def find_blast_files(self) -> List[Path]:
        """Find all .blast.tsv files in the input directory."""
        blast_files = list(self.input_dir.glob("*.blast.tsv"))
        if not blast_files:
            raise ValueError(f"No .blast.tsv files found in {self.input_dir}")
        
        logger.info(f"Found {len(blast_files)} BLAST files: {[f.name for f in blast_files]}")
        return blast_files
    
    def analyze_all_samples(self) -> List[pd.DataFrame]:
        """Analyze all BLAST files and return list of DataFrames."""
        blast_files = self.find_blast_files()
        dfs = []
        
        for blast_file in blast_files:
            logger.info(f"Analyzing {blast_file}")
            
            try:
                analyzer = BLASTCoverageAnalyzer(
                    blast_file=str(blast_file),
                    output_dir=self.output_dir,
                    top_n=self.top_n,
                )
                
                # Parse BLAST output
                df = analyzer.parse_blast_output()
                
                # Calculate coverage
                coverage_results = analyzer.calculate_reference_coverage(df)
                
                # Convert to DataFrame and add sample column
                df_sample = analyzer.create_summary_dataframe(coverage_results)
                dfs.append(df_sample)
                
            except Exception as e:
                logger.error(f"Failed to analyze {blast_file}: {e}")
                continue
        
        return dfs
    
    def create_aggregated_summary(self, all_dataframes: List[pd.DataFrame]) -> pd.DataFrame:
        """Create aggregated summary across all samples using pd.concat."""
        if not all_dataframes:
            return pd.DataFrame()
        
        # Concatenate all DataFrames
        aggregated_df = pd.concat(all_dataframes, ignore_index=True)
        
        # Sort by query, then by reference
        aggregated_df = aggregated_df.sort_values(['sample', 'reference'])
        
        return aggregated_df
    
    def write_aggregated_outputs(self, df_aggregated: pd.DataFrame) -> None:
        """Write aggregated output files."""
        logger.info("Writing aggregated output files...")
        
        # Create aggregated summary using pd.concat
        df_out = df_aggregated.rename(columns=dict(OUTPUT_COLUMNS))
        
        if df_out.empty:
            logger.warning("No data to write - all samples failed analysis")
            return
        
        # Write aggregated summary
        summary_file = self.output_dir / "blast_coverage_summary.tsv"
        df_out.to_csv(summary_file, sep="\t", index=False)
        logger.info(f"Aggregated summary written to {summary_file}")

    def write_multiqc_table(self, df_aggregated: pd.DataFrame) -> None:
        """Write MultiQC table.
        """
        table_comment = dedent("""
        # plot_type: 'table'
        # section_name: 'BLAST Coverage'
        # section_href: 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews'
        # description: 'Coverage and identity statistics from BLAST alignments.'
        # pconfig:
        #     namespace: 'BLAST'
        # headers:
        #     sample:
        #         title: 'Sample'
        #         description: 'Sample name'
        #         format: '{}'
        #     reference:
        #         title: 'Reference'
        #         description: 'Reference sequence ID'
        #         scale: False
        #         format: '{}'
        #     reference_title:
        #         title: 'Reference Title'
        #         description: 'Reference sequence title'
        #         scale: False
        #         format: '{}'
        #     avg_identity_covered:
        #         title: 'Avg Identity'
        #         description: 'Average BLAST percent identity. This is the average percent identity for all alignments to the same reference. The identity is weighted by the alignment length for each alignment.'
        #         format: '{:.1%}'
        #     reference_coverage:
        #         title: 'Reference Coverage'
        #         description: 'Reference coverage (proportion of reference covered by alignments)'
        #         format: '{:.1%}'
        #     total_alignment_length:
        #         title: 'Total Alignment Length'
        #         description: 'Total alignment length. This is the sum of the alignment lengths for all alignments to the same reference.'
        #         format: '{:.0f}'
        #     reference_length:
        #         title: 'Reference Length'
        #         description: 'Reference sequence length'
        #         format: '{:.0f}'
        #     num_alignments:
        #         title: '# Alignments'
        #         description: 'Number of BLAST alignments'
        #         format: '{:.0f}'
        """)
        multiqc_table_file = self.output_dir / "blast_coverage_mqc.txt"
        logger.info(f"Writing MultiQC table to {multiqc_table_file}")
        with open(multiqc_table_file, 'w') as fh:
            fh.write(table_comment)
            df_aggregated.to_csv(fh, sep="\t", index=False)
        logger.info(f"Wrote MultiQC table to '{multiqc_table_file}'")


def version_callback(value: bool):
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()


@app.command()
def main(
    input_dir: Path = typer.Argument(
        ...,
        help="Input directory containing .blast.tsv files (one per sample)"
    ),
    output_dir: Path = typer.Option(
        ".",
        "-o", "--output-dir",
        help="Output directory for results"
    ),
    top_n: int = typer.Option(
        1,
        "-n", "--top-n",
        help="Number of top references to analyze per query (default: 1)"
    ),
    verbose: bool = typer.Option(
        False,
        "-v", "--verbose",
        help="Enable verbose logging"
    ),
    version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
):
    """Analyze BLAST tabular output for coverage and identity statistics across multiple samples.
    
    The input directory should contain one or more .blast.tsv files, one for each sample.
    Reports are aggregations of top N results for each sample/tsv file.
    
    Examples:
      blast_coverage.py /path/to/blast/results/ -o results/
      blast_coverage.py /path/to/blast/results/ --top-n 3 -o analysis/
      blast_coverage.py /path/to/blast/results/ --verbose
    """
    init_logging(verbose)
    
    analyzer = MultiSampleBLASTAnalyzer(
        input_dir=input_dir,
        output_dir=output_dir,
        top_n=top_n,
    )
    
    all_dataframes = analyzer.analyze_all_samples()
    df_aggregated = analyzer.create_aggregated_summary(all_dataframes)
    analyzer.write_aggregated_outputs(df_aggregated)
    analyzer.write_multiqc_table(df_aggregated)
    logger.info("Done!")


if __name__ == "__main__":
    app()

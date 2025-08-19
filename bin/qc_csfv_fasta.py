#!/usr/bin/env python3
"""
CSFV FASTA Quality Analysis Tool

This script analyzes CSFV consensus FASTA files to provide comprehensive quality metrics
including N character distribution, GC content, reading frame analysis, and assembly quality.
"""

import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from collections import defaultdict
import statistics
from textwrap import dedent

import typer
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from rich.traceback import install
from Bio import SeqIO
from Bio.Seq import Seq
import edlib

# Version information
VERSION = "1.0.0"

# Install rich traceback handler for better error reporting
install(show_locals=True)

# Setup rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)

app = typer.Typer(help="Analyze CSFV FASTA consensus files for quality metrics")


def set_log_level(verbose: int):
    """Set logging level based on verbose flag."""
    if verbose == 0:
        logging.getLogger().setLevel(logging.INFO)
    elif verbose == 1:
        logging.getLogger().setLevel(logging.DEBUG)
    elif verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
        # Also enable debug for edlib and other libraries
        logging.getLogger("edlib").setLevel(logging.DEBUG)


# Polyprotein translation from NC_038912.1
CSFV_ALFORT_POLYPROTEIN = "MELNHFELLYKTNKQKPMGVEEPVYDATGKPLFGDPSEVHPQSTLKLPHDRGRGNIKTTLKNLPRKGDCRSGNHLGPVSGIYVKPGPVFYQDYMGPVYHRAPLEFFSEAQFCEVTKRIGRVTGSDGRLYHIYVCIDGCILLKLAKRGEPRTLKWIRNFTDCPLWVTSCSDDGASGSKEKKPDRINKGKLKIAPKEHEKDSRTKPPDATIVVEGVKYQVKKKGKVKGKNTQDGLYHNKNKPPESRKKLEKALLAWAVIAIMLYQPVEAENITQWNLSDNGTNGIQHAMYLRGISRSLHGIWPEKICKGVPTYLATDTELKEIQGMMDASEGTNYTCCKLQRHEWNKHGWCNWYNIDPWIQLMNRTQANLAEGPPAKECAVTCRYDKDADVNVVTQARNRPTTLTGCKKGKNFSFAGTVIEGPCNFNVSVEDILYGDHECGSLLQDTALYLVDGMTNTIENARQGAARVTSWLGRQLRIAGKRLEGRSKTWFGAYALSPYCNVTSKIGYIWYTNNCTPACLPKNTKIIGPGKFDTNAEDGKILHEMGGHLSEFLLLSLVVLSDFAPETASALYLILHYVIPQSHEEPEGCDTNQLNLTVELRTEDVIPSSVWNVGKYVCVRPDWWPYETKVALLFEEAGQVVKLALRALRDLTRVWNSASTTAFLICLIKVLRGQIVQGVIWLLLVTGAQGQLACKEDYRYAISSTNEIGLLGAGGLTTTWKEYNHDLQLNDGTVKAICVAGSFKVTALNVVSRRYLASLHKEALPTSVTFELLFDGTNPSTEEMGDDFGFGLCPFDTSPVVKGKYNTTLLNGSAFYLVCPIGWTGVIECTAVSPTTLRTEVVKTFRRDKPFPHRMDCVTTTVENEDLFYCKLGGNWTCVKGEPVVYTGGLVKQCRWCGFDFNEPDGLPHYPIGKCILANETGYRIVDSTDCNRDGVVISTEGSHECLIGNTTVKVHASDERLGPMPCRPKEIVSSAGPVRKTSCTFNYAKTLKNKYYEPRDSYFQQYMLKGEYQYWFDLDVTDRHSDYFAEFVVLVVVALLGGRYILWLIVTYIVLTEQLAAGLPLGQGEVVLIGNLITHTDIEVVVYFLLLYLVMRDEPIKKWILLLFHAMTNNPVKTITVALLMVSGVARGGKIDGGWQRLPETSFDIQLALTVIVVAVMLLAKRDPTTVPLVVTVATLRTAKMTNGLSTDIAIATVSTALLTWTYISDYYRYKTWLQYLISTVTGIFLIRVLKGIGELDLHTPTLPSYRPLFFILVYLISTAVVTRWNLDIAGLLLQCVPTLLMVFTMWADILTLILILPTYELTKLYYLKEVKIGAERGWLWKTNFKRVNDIYEVDQAGEGVYLFPSKQKTSSITGTMLPLIKAILISCISNKWQFIYLLYLIFEVSYYLHKKIIDEIAGGTNFISRLVAALIEANWAFDNEEVRGLKKFFLLSSRVKELIIKHKVRNEVMVHWFGDEEVYGMPKLVGLVKAATLSKNKHCILCTVCEDREWRGETCPKCGRFGPPMTCGMTLADFEEKHYKRIFFREDQSEGPVREEYAGYLQYRARGQLFLRNLPVLATKVKMLLVGNLGTEVGDLEHLGWVLRGPAVCKKVTEHEKCTTSIMDKLTAFFGVMPRGTTPRAPVRFPTSLLKIRRGLETGWAYTHQGGISSVDHVTCGKDLLVCDTMGRTRVVCQSNNKMTDESEYGVKTDSGCPEGARCYVFNPEAVNISGTKGAMVHLQKTGGEFTCVTASGTPAFFDLKNLKGWSGLPIFEASSGRVVGRVKVGKNEDSKPTKLMSGIQTVSKSTTDLTEMVKKITTMNRGEFRQITLATGAGKTTELPRSVIEEIGRHKRVLVLIPLRAAAESVYQYMRQKHPSIAFNLRIGEMKEGDMATGITYASYGYFCQMPQPKLRAAMVEYSFIFLDEYHCATPEQLAIMGKIHRFSENLRVVAMTATPAGTVTTTGQKHPIEEFIAPEVMKGEDLGSEYLDIAGLKIPVEEMKSNMLVFVPTRNMAVETAKKLKAKGYNSGYYYSGEDPSNLRVVTSQSPYVVVATNAIESGVTLPDLDVVVDTGLKCEKRIRLSPKMPFIVTGLKRMAVTIGEQAQRRGRVGRVKPGRYYRSQETPVGSKDYHYDLLQAQRYGIEDGINITKSFREMNYDWSLYEEDSLMITQLEILNNLLISEELPMAVKNIMARTDHPEPIQLAYNSYETQVPVLFPKIKNGEVTDSYDNYTFLNARKLGDDVPPYVYATEDEDLAVELLGLDWPDPGNQGTVEAGRALKQVVGLSTAENALLVALFGYVGYQALSKRHIPVVTDIYSIEDHRLEDTTHLQYAPNAIKTEGKETELKELAQGDVQRCVEAMTNYAREGIQFMKSQALKVKETPTYKETMDTVTDYVKKFMEALTDSKEDIIKYGLWGTHTALYKSICARLGSETAFATLVVKWLAFGGESIADHVKQAATDLVVYYIINRPQFPGDTETQQEGRKFVASLLVSALVTYTYKSWNYNNLSKIVEPALATLPYAATALKLFAPTRLESVVILSTAIYKTYLSIRRGKSDGLLGTGVSAAMEIMSQNPVSVGIAVMLGVGAVAAHNAIEASEQKRTLLMKVFVKNFLDQAATDELVKESPEKIIMALFEAVQTVGNPLRLVYHLYGVFYKGWEAKELAQRTAGRNLFTLIMFEAVELLGVDSEGKIRQLSSNYILELLYKFRDSIKSSVREMAISWAPAPFSCDWTPTDDRIGLPQDNFLQVETKCPCGYKMKAVKNCAGELRLLEEEGSFLCRNKFGRGSRNYRVTKYYDDNLSEIKPVIRMEGHVELYYKGATIKLDFNNSKTILATDKWEVDHSTLVRVLKRHTGAGYHGAYLGEKPNHKHLIERDCATITKDKVCFLKMKRGCAFTYDLSLHNLTRLIELVHKNNLEDKEIPAVTVTTWLAYTFVNEDIGTIKPAFGEKVTPEMQEEITLQPAVVVDTTDVTVTVVGEAPTMTTGETPTAFTSSGSDPKGQQVLKLGVGEGQYPGTNPQRASLHEAIQGADERPSVLILGSDKATSNRVKTAKNVKVYRGRDPLEVRDMMRRGKILVIALSRVDNALLKFVDYKGTFLTRETLEALSLGRPKKKNITKAEAQWLLCLEDQMEELPDWFAAGEPIFLEANIKHDRYHLVGDIATIKEKAKQLGATDSTKISKEVGAKVYSMKLSNWVMQEENKQGNLTPLFEELLQQCPPGGQNKTAHMVSAYQLAQGNWMPTSCHVFMGTISARRTKTHPYEAYVKLRELVEEHKMKTLCPGSSLGKHNEWIIGKIKYQGNLRTKHMLNPGKVAEQLCREGHRRNVYNKTIGSVMTATGIRLEKLPVVRAQTDTTNFHQAIRDKIDKEENLQTPGLHKKLMEVFNALKRPELESSYDAVEWEELERGINRKGAAGFFERKNIGEILDSEKNKVEEIIDNLKKGRNIKYYETAIPKNEKRDVNDDWTSGDFVDEKKPRVIQYPEAKTRLAITKVMYKWVKQKPVVIPGYEGKTPLFQIFDKVKKEWDQFQNPVAVSFDTKAWDTQVTTKDLELIKDIQKYYFKKKWHKFIDTLTMHMSEVPVISADGEVYIRKGQRGSGQPDTSAGNSMLNVLTMVYAFCEATGVPYKSFDRVAKIHVCGDDGFLITERALGEKFASKGVQILYEAGKPQKITEGDKMKVAYQFDDIEFCSHTPIQVRWSDNTSSYMPGRNTTTILAKMATRLDSSGERGTIAYEKAVAFSFLLMYSWNPLIRRICLLVLSTELQVKPGKSTTYYYEGDPISAYKEVIGHNLFDLKRTSFEKLAKLNLSMSVLGAWTRHTSKRLLQDCVNMGVKEGNWLVNADRLVSSKTGNRYIPGEGHTLQGRHYEELVLARKQINNFQGTDRYNLGPIVNMVLRRLRVMMMTLIGRGV"

COLUMNS = [
    ("sequence_id", "Sequence ID"),
    ("quality_comment", "QC Comment"),
    ("frameshift_evidence", "Frameshift Evidence"),
    ("coverage_percentage", "Genome Coverage (%)"),
    ("length", "Sequence Length"),
    ("total_ns", "Total Ns"),
    ("end_5prime_ns", "End 5' Ns"),
    ("end_3prime_ns", "End 3' Ns"),
    ("n_regions", "N Regions"),
    ("gc_content", "GC Content (%)"),
    ("ambiguous_bases", "Ambiguous Bases"),
    ("ambiguous_detail", "Ambiguous Detail"),
    ("polyprotein_start", "Polyprotein Start"),
    ("polyprotein_end", "Polyprotein End"),
    ("potential_frameshift", "Potential Frameshift"),
    ("stop_codons", "Stop Codons"),
    ("stop_positions", "Stop Positions"),
    ("best_alignment_edit_distance", "Best Alignment Edit Distance"),
    ("best_alignment_score", "Best Alignment Score"),
    ("homopolymer_issues", "Homopolymer Issues"),
]


def find_n_regions(sequence: str) -> List[Tuple[int, int, int]]:
    """Find all regions of N characters in the sequence."""
    n_regions = []
    in_n_region = False
    start_pos: Optional[int] = None

    for i, base in enumerate(sequence.upper()):
        if base == "N":
            if not in_n_region:
                start_pos = i + 1  # 1-based coordinate
                in_n_region = True
        else:
            if in_n_region and start_pos is not None:  # Check for None
                end_pos = i  # 1-based coordinate (end position of last N)
                length = end_pos - start_pos + 1
                n_regions.append((start_pos, end_pos, length))
                in_n_region = False

    # Handle case where sequence ends with Ns
    if in_n_region and start_pos is not None:  # Check for None
        end_pos = len(sequence)
        length = end_pos - start_pos + 1
        n_regions.append((start_pos, end_pos, length))

    return n_regions


def format_n_regions(n_regions: List[Tuple[int, int, int]]) -> str:
    """Format N regions as start..end strings separated by semicolons."""
    if not n_regions:
        return ""

    formatted = []
    for start, end, length in n_regions:
        if length == 1:
            formatted.append(str(start))
        else:
            formatted.append(f"{start}..{end}")

    return ";".join(formatted)


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content, ignoring N and ambiguous bases."""
    clean_seq = sequence.upper()
    valid_bases = sum(1 for base in clean_seq if base in "ATGC")

    if valid_bases == 0:
        return 0.0

    gc_count = sum(1 for base in clean_seq if base in "GC")
    return (gc_count / valid_bases) * 100


def count_ambiguous_bases(sequence: str) -> Dict[str, int]:
    """Count ambiguous bases (not ATGCN)."""
    ambiguous_bases = defaultdict(int)
    standard_bases = set("ATGCN")

    for base in sequence.upper():
        if base not in standard_bases:
            ambiguous_bases[base] += 1

    return dict(ambiguous_bases)


def find_polyprotein_boundaries(
    sequence: str, reference_polyprotein: str
) -> Tuple[int, int, str]:
    """Find the start and end of the polyprotein using edlib alignment with reference."""
    # Convert sequence to BioPython Seq object for translation
    seq_obj = Seq(sequence)

    logger.debug("=== Polyprotein Boundary Detection ===")
    logger.debug(f"Sequence length: {len(sequence)}")
    logger.debug(f"Reference polyprotein length: {len(reference_polyprotein)}")

    # Try all 6 reading frames and find the best alignment
    best_alignment = None
    best_score = -1
    best_protein = ""

    # For CSFV sequences, only consider forward strand (reverse=False)
    # This prevents false reverse strand detection
    for reverse in [False]:  # Only forward strand
        strand_name = "forward"  # Always forward for CSFV
        for frame in range(3):
            if reverse:
                # Reverse complement
                current_seq = seq_obj.reverse_complement()
            else:
                current_seq = seq_obj

            # Translate from the specified frame
            try:
                protein = current_seq[frame:].translate(
                    table="1", cds=False
                )  # Standard genetic code
                protein_str = str(protein)

                # Use edlib to align against reference polyprotein
                result = edlib.align(
                    protein_str, reference_polyprotein, mode="HW", task="path"
                )

                if result["editDistance"] != -1:  # Valid alignment
                    # Calculate alignment score (lower edit distance = better)
                    base_score = len(reference_polyprotein) - result["editDistance"]

                    # For CSFV, prefer forward strand unless reverse is significantly better
                    if reverse:
                        # Only consider reverse if it's significantly better than any forward alignment
                        # This prevents false reverse strand detection for CSFV
                        if base_score > 100:  # Only if alignment is reasonably good
                            score = base_score - 5000  # Still penalize reverse
                        else:
                            score = -10000  # Reject poor reverse alignments
                    else:
                        score = base_score

                    logger.debug(
                        f"  Frame {frame} ({strand_name}): Protein length {len(protein_str)}, "
                        f"Edit distance: {result['editDistance']}, Base score: {base_score}, "
                        f"Final score: {score}"
                    )

                    if score > best_score:
                        best_score = score
                        best_alignment = result
                        best_protein = protein_str

            except Exception:
                continue

    if best_alignment is None:
        return 0, 0, ""

    # Get the alignment locations from the best alignment
    locations = best_alignment["locations"]
    if not locations:
        return 0, 0, best_protein

    # Find the best polyprotein start by testing ATG positions in each reading frame
    # This gives us biologically accurate boundaries instead of hardcoded values

    best_start_pos = 0
    best_start_score = -1

    # Test each reading frame
    for frame in range(3):
        # Find all ATG positions in this frame, handling ambiguous bases
        # Ambiguous bases that could form ATG: A, T, G, M (A/C), R (A/G), W (A/T), Y (C/T)
        atg_positions = []
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i : i + 3].upper()
            # Check if this could be ATG, accounting for ambiguous bases
            if (
                codon[0] in "AMRW"  # A or ambiguous bases that could be A
                and codon[1] in "TYW"  # T or ambiguous bases that could be T
                and codon[2] in "GKR"
            ):  # G or ambiguous bases that could be G
                atg_positions.append(i)

        # Test each ATG position in this frame
        for atg_pos in atg_positions:
            # Extract sequence from this ATG onwards
            if atg_pos + len(reference_polyprotein) * 3 <= len(sequence):
                test_seq = sequence[atg_pos:]
                test_protein = Seq(test_seq).translate(table="1", cds=False)

                # Align this protein against reference
                result = edlib.align(
                    str(test_protein), reference_polyprotein, mode="HW", task="path"
                )
                if result["editDistance"] != -1:
                    score = len(reference_polyprotein) - result["editDistance"]
                    if score > best_start_score:
                        best_start_score = score
                        best_start_pos = atg_pos
    best_fallback_pos = 0
    best_fallback_score = float("inf")
    # If we found a good start position, use it
    if best_start_pos > 0:
        polyprotein_start = best_start_pos + 1  # Convert to 1-based
        polyprotein_end = polyprotein_start + len(reference_polyprotein) * 3 - 1

        # Ensure we don't go beyond sequence length
        if polyprotein_end > len(sequence):
            polyprotein_end = len(sequence)
            polyprotein_start = polyprotein_end - len(reference_polyprotein) * 3 + 1
    else:
        # Fallback: systematically test all possible start positions
        # Find the position with minimum edlib edit distance - completely data-driven
        # Test every position that could potentially start a polyprotein
        # Use a reasonable step size to avoid testing every single position
        step_size = 3  # Test every 3rd position (codon-aligned)

        for pos in range(0, len(sequence) - len(reference_polyprotein) * 3, step_size):
            test_seq = sequence[pos:]
            test_protein = Seq(test_seq).translate(table="1", cds=False)

            # Align against reference
            result = edlib.align(
                str(test_protein), reference_polyprotein, mode="HW", task="path"
            )
            if (
                result["editDistance"] != -1
                and result["editDistance"] < best_fallback_score
            ):
                best_fallback_score = result["editDistance"]
                best_fallback_pos = pos

        # If we found a good position, use it
        if best_fallback_pos > 0:
            polyprotein_start = best_fallback_pos + 1  # Convert to 1-based
            polyprotein_end = polyprotein_start + len(reference_polyprotein) * 3 - 1
            if polyprotein_end > len(sequence):
                polyprotein_end = len(sequence)
                polyprotein_start = polyprotein_end - len(reference_polyprotein) * 3 + 1
        else:
            # Last resort: use sequence length to determine boundaries
            polyprotein_start = 1
            polyprotein_end = len(sequence)

    if best_start_pos > 0:
        logger.debug(
            f"  Found best polyprotein start at ATG position {best_start_pos} (score: {best_start_score})"
        )
        logger.debug(
            f"  Dynamic polyprotein boundaries: {polyprotein_start}..{polyprotein_end}"
        )
    else:
        if best_fallback_pos > 0:
            logger.debug(
                f"  Fallback: systematic search found best position {best_fallback_pos} (edit distance: {best_fallback_score})"
            )
            logger.debug(
                f"  Fallback polyprotein boundaries: {polyprotein_start}..{polyprotein_end}"
            )
        else:
            logger.debug(
                "  Fallback: no good position found, using full sequence boundaries"
            )
            logger.debug(
                f"  Fallback polyprotein boundaries: {polyprotein_start}..{polyprotein_end}"
            )
    logger.debug(f"  Best alignment score: {best_score}")
    logger.debug("=" * 80)

    return polyprotein_start, polyprotein_end, best_protein


def score_reading_frame_with_polyprotein(
    protein: str, polyprotein_start: int, polyprotein_end: int
) -> float:
    """Score a reading frame considering only the polyprotein region."""
    if not protein:
        return float("inf")

    # Only consider stop codons within the polyprotein region
    stop_positions = [i for i, aa in enumerate(protein) if aa == "*"]

    if not stop_positions:
        return 0.0  # Perfect - no stops

    # Filter stops to only those within the polyprotein region
    polyprotein_stops = [
        pos for pos in stop_positions if polyprotein_start <= pos < polyprotein_end
    ]

    if not polyprotein_stops:
        return 0.0  # No stops in polyprotein region

    # Score based on position within polyprotein (early stops are worse)
    early_stops = [pos for pos in polyprotein_stops if pos < polyprotein_start + 100]
    mid_stops = [
        pos
        for pos in polyprotein_stops
        if polyprotein_start + 100 <= pos < polyprotein_start + 1000
    ]
    late_stops = [pos for pos in polyprotein_stops if pos >= polyprotein_start + 1000]

    score = len(early_stops) * 10.0 + len(mid_stops) * 2.0 + len(late_stops) * 0.5

    return score


def find_best_reading_frame_with_polyprotein(
    sequence: str, reference_polyprotein: str
) -> Tuple[int, str, List[int], bool, int, int]:
    """Find the best reading frame considering polyprotein boundaries."""
    # First find polyprotein boundaries
    poly_start, poly_end, best_protein = find_polyprotein_boundaries(
        sequence, reference_polyprotein
    )

    if not best_protein:
        return 0, "", [], False, 0, 0

    # Determine which frame and strand gave us the best alignment
    best_frame = 0
    best_strand = False
    best_score = float("inf")
    best_stop_positions = []

    logger.debug(f"=== Reading Frame Analysis for sequence length {len(sequence)} ===")
    logger.debug(f"Polyprotein boundaries: {poly_start}..{poly_end}")

    # Check all 6 reading frames (only forward strand for CSFV)
    for reverse in [False]:  # Only forward strand for CSFV
        strand_name = "forward"  # Always forward for CSFV
        for frame in range(3):
            # Inline the translation logic instead of calling translate_sequence_biopython
            seq_obj = Seq(sequence)
            if reverse:
                seq_obj = seq_obj.reverse_complement()
            try:
                protein = seq_obj[frame:].translate(table="1", cds=False)
                protein = str(protein)
            except Exception:
                protein = ""

            if protein:
                # Calculate edlib edit distance BEFORE trimming
                full_alignment = edlib.align(
                    protein, reference_polyprotein, mode="HW", task="path"
                )
                full_edit_distance = full_alignment["editDistance"]

                # Calculate edlib edit distance AFTER trimming to polyprotein region
                poly_protein_start = poly_start - 1  # Convert to 0-based
                poly_protein_end = poly_end
                trimmed_protein = protein[poly_protein_start:poly_protein_end]

                trimmed_alignment = edlib.align(
                    trimmed_protein, reference_polyprotein, mode="HW", task="path"
                )
                trimmed_edit_distance = trimmed_alignment["editDistance"]

                logger.debug(
                    f"Frame {frame} ({strand_name}): Full protein length: {len(protein)}, "
                    f"Full edit distance: {full_edit_distance}, "
                    f"Trimmed protein length: {len(trimmed_protein)}, "
                    f"Trimmed edit distance: {trimmed_edit_distance}"
                )

                score = score_reading_frame_with_polyprotein(
                    protein, poly_start, poly_end
                )

                if score < best_score:
                    best_score = score
                    best_frame = frame
                    best_strand = reverse
                    best_protein = protein
                    stop_positions = [i for i, aa in enumerate(protein) if aa == "*"]
                    best_stop_positions = stop_positions

    logger.debug(
        f"Selected: Frame {best_frame} ({'reverse' if best_strand else 'forward'}) with score {best_score}"
    )
    logger.debug("=" * 80)

    return (
        best_frame,
        best_protein,
        best_stop_positions,
        best_strand,
        poly_start,
        poly_end,
    )


def find_homopolymer_regions(
    sequence: str, min_length: int = 5
) -> List[Tuple[int, int, int, str]]:
    """Find homopolymer regions (stretches of the same base) in a DNA sequence.

    Args:
        sequence: DNA sequence to analyze
        min_length: Minimum length to consider a homopolymer (default: 5)

    Returns:
        List of tuples: (start_pos, end_pos, length, base_type)
        Positions are 1-based for consistency with other functions
    """
    homopolymers = []
    current_base = None
    current_start = 0
    current_length = 0

    for i, base in enumerate(sequence.upper()):
        if base in "ATGC":
            if base == current_base:
                current_length += 1
            else:
                # Check if previous homopolymer meets minimum length
                if current_length >= min_length:
                    homopolymers.append(
                        (
                            current_start + 1,  # Convert to 1-based
                            current_start + current_length,  # Convert to 1-based
                            current_length,
                            current_base,
                        )
                    )

                # Start new homopolymer
                current_base = base
                current_start = i
                current_length = 1
        else:
            # Non-standard base, reset homopolymer
            if current_length >= min_length:
                homopolymers.append(
                    (
                        current_start + 1,  # Convert to 1-based
                        current_start + current_length,  # Convert to 1-based
                        current_length,
                        current_base,
                    )
                )

            current_base = None
            current_length = 0

    # Check final homopolymer
    if current_length >= min_length:
        homopolymers.append(
            (
                current_start + 1,  # Convert to 1-based
                current_start + current_length,  # Convert to 1-based
                current_length,
                current_base,
            )
        )

    return homopolymers


def analyze_homopolymers_and_stops(
    sequence: str, stop_positions: List[int], poly_start: int, poly_end: int
) -> Dict[str, Any]:
    """Analyze homopolymer regions in relation to premature stop codons.

    This function systematically checks for homopolymers upstream of ALL premature
    stop codons to identify potential sequencing errors causing frameshifts.

    Args:
        sequence: DNA sequence to analyze
        stop_positions: List of stop codon positions (1-based)
        poly_start: Start of polyprotein region (1-based)
        poly_end: End of polyprotein region (1-based)

    Returns:
        Dictionary with comprehensive analysis of homopolymer-stop relationships
    """
    if not stop_positions:
        return {
            "has_homopolymer_issues": False,
            "primary_issue": None,
            "homopolymer_details": None,
            "analysis_summary": "No premature stop codons detected",
        }

    # Find stop codons in the polyprotein region
    polyprotein_stops = [pos for pos in stop_positions if poly_start <= pos <= poly_end]

    if not polyprotein_stops:
        return {
            "has_homopolymer_issues": False,
            "primary_issue": None,
            "homopolymer_details": None,
            "analysis_summary": "No premature stop codons in polyprotein region",
        }

    # Sort stops by position
    sorted_stops = sorted(polyprotein_stops)
    first_stop = sorted_stops[0]

    # Check for early stops (within first 100bp of polyprotein)
    early_stops = [pos for pos in sorted_stops if pos < poly_start + 100]

    # Look for homopolymers upstream of each stop codon
    homopolymer_issues = []

    for stop_pos in sorted_stops:
        # Search up to 100 nucleotides upstream of each stop (frameshifts can propagate from further upstream)
        search_start = max(poly_start, stop_pos - 100)
        search_end = stop_pos

        # Extract the sequence region to search
        region_seq = sequence[search_start - 1 : search_end]  # Convert to 0-based

        # Find homopolymers in this region (3bp+ can cause frameshifts, but 5bp+ are more likely)
        homopolymers = find_homopolymer_regions(region_seq, min_length=3)

        for homopolymer in homopolymers:
            # Calculate distance from homopolymer end to stop
            homopolymer_end_in_full = (
                search_start + homopolymer[1] - 1
            )  # Convert back to 1-based
            distance_to_stop = stop_pos - homopolymer_end_in_full

            # Only consider homopolymers that are upstream of the stop
            if (
                0 <= distance_to_stop <= 50
            ):  # Within 50 nucleotides upstream (frameshifts can propagate)
                homopolymer_issues.append(
                    {
                        "start": homopolymer_end_in_full
                        - homopolymer[2]
                        + 1,  # Start position
                        "end": homopolymer_end_in_full,  # End position
                        "length": homopolymer[2],
                        "base": homopolymer[3],
                        "distance_to_stop": distance_to_stop,
                        "stop_position": stop_pos,
                        "is_early_stop": stop_pos in early_stops,
                    }
                )

    # Analyze the findings
    if not homopolymer_issues:
        # No homopolymers found, but we have stops - analyze why
        if early_stops:
            return {
                "has_homopolymer_issues": False,
                "primary_issue": "early_stops",
                "homopolymer_details": None,
                "analysis_summary": f"Early stop codons at positions {early_stops} - may indicate assembly issues or true biological variation",
            }
        else:
            return {
                "has_homopolymer_issues": False,
                "primary_issue": "distributed_stops",
                "homopolymer_details": None,
                "analysis_summary": f"Stop codons distributed throughout polyprotein (positions: {sorted_stops[:5]}{'...' if len(sorted_stops) > 5 else ''}) - may indicate frameshift without clear homopolymer cause",
            }

    # Check if there's a homopolymer upstream of the FIRST stop codon
    first_stop_homopolymers = [
        h for h in homopolymer_issues if h["stop_position"] == first_stop
    ]

    # Initialize most_significant
    most_significant = None

    if first_stop_homopolymers:
        # Found homopolymer upstream of first stop - this is the most important
        most_significant = max(
            first_stop_homopolymers, key=lambda x: (x["length"], -x["distance_to_stop"])
        )
        analysis_summary = (
            f"First premature stop codon at position {most_significant['stop_position']} "
            f"possibly due to sequencing errors in homopolymer region of "
            f"{most_significant['length']}bp {most_significant['base']} "
            f"at position {most_significant['start']}-{most_significant['end']} "
            f"({most_significant['distance_to_stop']}bp upstream)"
        )
    else:
        # No homopolymer upstream of first stop - check other stops for homopolymers
        if homopolymer_issues:
            # Found homopolymers upstream of other stops
            most_significant = max(
                homopolymer_issues,
                key=lambda x: (x["is_early_stop"], x["length"], -x["distance_to_stop"]),
            )

            if most_significant["is_early_stop"]:
                analysis_summary = (
                    f"First premature stop codon at position {first_stop} - no homopolymer region found upstream "
                    f"(searched 100bp upstream). However, early stop codon at position {most_significant['stop_position']} "
                    f"possibly due to sequencing errors in homopolymer region of "
                    f"{most_significant['length']}bp {most_significant['base']} "
                    f"at position {most_significant['start']}-{most_significant['end']} "
                    f"({most_significant['distance_to_stop']}bp upstream)"
                )
            else:
                analysis_summary = (
                    f"First premature stop codon at position {first_stop} - no homopolymer region found upstream "
                    f"(searched 100bp upstream). However, frameshift possibly due to sequencing errors in homopolymer region of "
                    f"{most_significant['length']}bp {most_significant['base']} "
                    f"at position {most_significant['start']}-{most_significant['end']} "
                    f"upstream of stop codon at position {most_significant['stop_position']} "
                    f"({most_significant['distance_to_stop']}bp upstream)"
                )
        else:
            # No homopolymers found at all
            if first_stop < poly_start + 100:  # First stop is early
                analysis_summary = (
                    f"First premature stop codon at position {first_stop} - no homopolymer region found upstream "
                    f"(searched 100bp upstream). This may indicate assembly issues or true biological variation."
                )
            else:
                analysis_summary = (
                    f"First premature stop codon at position {first_stop} - no homopolymer region found upstream "
                    f"(searched 100bp upstream)."
                )

    return {
        "has_homopolymer_issues": True,
        "primary_issue": "homopolymer_related",
        "homopolymer_details": most_significant,  # Now guaranteed to be defined
        "analysis_summary": analysis_summary,
    }


def count_stop_codons_in_polyprotein(
    sequence: str, frame: int, reverse: bool, poly_start: int, poly_end: int
) -> List[int]:
    """Count stop codons only within the polyprotein region."""
    if reverse:
        seq_obj = Seq(sequence).reverse_complement()
    else:
        seq_obj = Seq(sequence)

    # Translate the polyprotein region only
    poly_nuc_start = poly_start - 1  # Convert to 0-based
    poly_nuc_end = poly_end

    # Extract the polyprotein nucleotide sequence
    poly_nuc_seq = seq_obj[poly_nuc_start:poly_nuc_end]

    # Translate this region
    poly_protein = poly_nuc_seq.translate(table="1", cds=False)

    logger.debug(
        f"  Polyprotein region {poly_start}..{poly_end} (nuc: {poly_nuc_start}..{poly_nuc_end})"
    )
    logger.debug(f"  Translated protein length: {len(poly_protein)}")
    logger.debug(f"  First 50 AA: {str(poly_protein[:50])}")
    logger.debug(f"  Last 50 AA: {str(poly_protein[-50:])}")

    # Find stop codons in the polyprotein protein sequence
    stop_positions = []
    for i, aa in enumerate(poly_protein):
        if aa == "*":
            # Convert back to nucleotide coordinates in original sequence
            if reverse:
                nt_pos = len(sequence) - (poly_nuc_start + (i * 3))
            else:
                nt_pos = poly_nuc_start + (i * 3) + 1
            stop_positions.append(nt_pos)

    logger.debug(f"  Found {len(stop_positions)} stop codons in polyprotein region")
    if stop_positions:
        logger.debug(f"  First 10 stop positions: {stop_positions[:10]}")

    return stop_positions


def analyze_reading_frame(sequence: str) -> Dict[str, Any]:
    """Analyze reading frame and detect potential frameshifts using reference polyprotein alignment."""
    frame, protein, _, reverse, poly_start, poly_end = (
        find_best_reading_frame_with_polyprotein(sequence, CSFV_ALFORT_POLYPROTEIN)
    )

    if not protein:
        # Fallback if polyprotein alignment fails
        return {
            "best_frame": 1,
            "strand": "forward",
            "total_stops": 0,
            "stop_positions": [],
            "protein_length": 0,
            "start_codons": 0,
            "start_codon_positions": [],
            "polyprotein_start": 0,
            "polyprotein_end": 0,
            "potential_frameshift": False,
            "frameshift_evidence": "Polyprotein alignment failed",
        }

    # Now count stop codons ONLY within the polyprotein region
    polyprotein_stop_positions = count_stop_codons_in_polyprotein(
        sequence, frame, reverse, poly_start, poly_end
    )

    # Check for homopolymer regions that might cause sequencing errors
    homopolymer_analysis = analyze_homopolymers_and_stops(
        sequence, polyprotein_stop_positions, poly_start, poly_end
    )

    # Find start codons in the best reading frame
    # Inline the start codon detection logic instead of calling find_start_codons_biopython
    start_codons = []
    seq_obj = Seq(sequence)
    if reverse:
        seq_obj = seq_obj.reverse_complement()
    for i in range(frame, len(seq_obj) - 2, 3):
        codon = str(seq_obj[i : i + 3]).upper()
        # Check if this could be ATG, accounting for ambiguous bases
        if (
            codon[0] in "AMRW"  # A or ambiguous bases that could be A
            and codon[1] in "TYW"  # T or ambiguous bases that could be T
            and codon[2] in "GKR"
        ):  # G or ambiguous bases that could be G
            start_codons.append(i)

    start_codon_positions = []
    if reverse:
        for pos in start_codons:
            start_codon_positions.append(len(sequence) - pos)
    else:
        start_codon_positions = [pos + 1 for pos in start_codons]  # 1-based

    # Get alignment information for the chosen polyprotein start
    # Inline the alignment info logic instead of calling get_alignment_info
    alignment_info = {"edit_distance": -1, "score": -1}
    if poly_start > 0:
        # Extract sequence from polyprotein start onwards
        if poly_start + len(CSFV_ALFORT_POLYPROTEIN) * 3 <= len(sequence):
            test_seq = sequence[poly_start - 1 :]  # Convert to 0-based
            test_protein = Seq(test_seq).translate(table="1", cds=False)

            # Align against reference
            result = edlib.align(
                str(test_protein), CSFV_ALFORT_POLYPROTEIN, mode="HW", task="path"
            )
            if result["editDistance"] != -1:
                score = len(CSFV_ALFORT_POLYPROTEIN) - result["editDistance"]
                alignment_info = {
                    "edit_distance": result["editDistance"],
                    "score": score,
                }

    analysis = {
        "best_frame": frame + 1,  # 1-based for reporting
        "strand": "reverse" if reverse else "forward",
        "total_stops": len(polyprotein_stop_positions),  # Only stops in polyprotein
        "stop_positions": polyprotein_stop_positions,
        "protein_length": len(protein),
        "start_codons": len(start_codons),
        "start_codon_positions": start_codon_positions,
        "polyprotein_start": poly_start,
        "polyprotein_end": poly_end,
        "best_alignment_edit_distance": alignment_info["edit_distance"],
        "best_alignment_score": alignment_info["score"],
        "homopolymer_issues": str(homopolymer_analysis["homopolymer_details"])
        if homopolymer_analysis["homopolymer_details"]
        else "",
    }

    # Improved frameshift detection considering polyprotein boundaries
    if len(polyprotein_stop_positions) > 0:
        # Use the comprehensive homopolymer analysis
        if homopolymer_analysis["has_homopolymer_issues"]:
            analysis["potential_frameshift"] = True
            analysis["frameshift_evidence"] = homopolymer_analysis["analysis_summary"]
            return analysis
        else:
            # No homopolymer issues found, but we still have stops
            analysis["potential_frameshift"] = True
            analysis["frameshift_evidence"] = homopolymer_analysis["analysis_summary"]
            return analysis
    else:
        analysis["potential_frameshift"] = False
        analysis["frameshift_evidence"] = "No stop codons detected"

    return analysis


def generate_quality_comment(
    sequence: str,
    n_regions: List[Tuple[int, int, int]],
    ambiguous_bases: Dict[str, int],
    reading_frame_analysis: Dict,
    total_ns: int,
) -> str:
    """Generate comprehensive quality assessment comment."""
    seq_length = len(sequence)
    comments = []

    # Genome coverage comment
    coverage_percentage = (
        ((seq_length - total_ns) / seq_length * 100) if seq_length > 0 else 0
    )
    if coverage_percentage > 99:
        coverage_comment = (
            f"High reference genome coverage of {coverage_percentage:.1f}%"
        )
    elif coverage_percentage > 90:
        coverage_comment = (
            f"Moderate reference genome coverage of {coverage_percentage:.1f}%"
        )
    else:
        coverage_comment = f"Low reference genome coverage of {coverage_percentage:.1f}% possibly due to sequencing, sample or PCR amplification issues"
    comments.append(coverage_comment)

    n_comment = []
    # Total Ns
    if total_ns > 0:
        n_comment.append(f"{total_ns} Ns overall")
    else:
        n_comment.append("No Ns detected")

    # N regions at ends
    if n_regions:
        first_region = n_regions[0]
        last_region = n_regions[-1]

        # 5' end Ns
        if first_region[0] == 1:
            length = first_region[2]
            n_comment.append(f"{length}bp Ns at 5' end")

        # 3' end Ns
        if last_region[1] == seq_length:
            length = last_region[2]
            n_comment.append(f"{length}bp Ns at 3' end")

    # Ambiguous bases
    if ambiguous_bases:
        total_ambiguous = sum(ambiguous_bases.values())
        n_comment.append(f"{total_ambiguous} ambiguous bases")
    comments.append(". ".join(n_comment))

    # Premature stop codons (first position only, no numbers)
    if reading_frame_analysis["total_stops"] > 0:
        stop_positions = reading_frame_analysis.get("stop_positions", [])
        if stop_positions:
            first_stop = min(stop_positions)
            comments.append(f"First premature stop codon detected at {first_stop}")

            # Add reason for stop codon if available
            frameshift_evidence = reading_frame_analysis.get("frameshift_evidence", "")
            if frameshift_evidence and "homopolymer" in frameshift_evidence.lower():
                hp_comment_index = frameshift_evidence.find("homopolymer region")
                homopolymer_region = frameshift_evidence[hp_comment_index:]
                comments.append(
                    f"Possibly due to sequencing errors in {homopolymer_region}"
                )
            elif frameshift_evidence and "early stop" in frameshift_evidence.lower():
                comments.append("Early stop codon may indicate assembly issues")
    else:
        comments.append("No premature stop codons detected")

    # Generate final comment
    if not comments:
        return "High quality assembly with little to no issues"

    comment = ". ".join(comments) + "."
    return comment


def analyze_fasta_file(filepath: Path) -> List[Dict[str, Any]]:
    """Analyze a single FASTA file and return summary statistics."""
    results = []

    try:
        for record in SeqIO.parse(filepath, "fasta"):
            seq_id = record.id
            sequence = str(record.seq)

            # Basic metrics
            seq_length = len(sequence)
            total_ns = sequence.upper().count("N")

            # N regions analysis
            n_regions = find_n_regions(sequence)
            n_regions_str = format_n_regions(n_regions)

            # End N analysis
            end_ns = {"5_prime": 0, "3_prime": 0}
            if n_regions:
                if n_regions[0][0] == 1:  # Starts with Ns
                    end_ns["5_prime"] = n_regions[0][2]
                if n_regions[-1][1] == seq_length:  # Ends with Ns
                    end_ns["3_prime"] = n_regions[-1][2]

            # GC content and ambiguous bases
            gc_content = calculate_gc_content(sequence)
            ambiguous_bases = count_ambiguous_bases(sequence)

            # Reading frame analysis
            reading_frame_analysis = analyze_reading_frame(sequence)

            # Coverage metrics
            coverage_percentage = (
                ((seq_length - total_ns) / seq_length * 100) if seq_length > 0 else 0
            )
            qc_comment = generate_quality_comment(
                sequence, n_regions, ambiguous_bases, reading_frame_analysis, total_ns
            )
            result = {
                "sequence_id": seq_id,
                "length": seq_length,
                "total_ns": total_ns,
                "coverage_percentage": coverage_percentage,
                "n_regions": n_regions_str,
                "end_5prime_ns": end_ns["5_prime"],
                "end_3prime_ns": end_ns["3_prime"],
                "gc_content": gc_content,
                "ambiguous_bases": sum(ambiguous_bases.values()),
                "ambiguous_detail": ambiguous_bases,
                "best_reading_frame": reading_frame_analysis["best_frame"],
                "strand": reading_frame_analysis["strand"],
                "start_codons": reading_frame_analysis["start_codons"],
                "stop_codons": reading_frame_analysis["total_stops"],
                "potential_frameshift": reading_frame_analysis["potential_frameshift"],
                "frameshift_evidence": reading_frame_analysis.get(
                    "frameshift_evidence", ""
                ),
                "polyprotein_start": reading_frame_analysis["polyprotein_start"],
                "polyprotein_end": reading_frame_analysis["polyprotein_end"],
                "best_alignment_edit_distance": reading_frame_analysis.get(
                    "best_alignment_edit_distance", -1
                ),
                "best_alignment_score": reading_frame_analysis.get(
                    "best_alignment_score", -1
                ),
                "stop_positions": ";".join(
                    map(str, reading_frame_analysis.get("stop_positions", []))
                ),
                "homopolymer_issues": reading_frame_analysis.get(
                    "homopolymer_issues", ""
                ),
                "quality_comment": qc_comment,
            }

            results.append(result)

    except Exception as e:
        logger.error(f"Error analyzing {filepath}: {e}")
        raise

    return results


def write_multiqc_table(df: pd.DataFrame, output_file: Path) -> None:
    """Write MultiQC table with YAML comment header."""
    
    # Drop columns that shouldn't be in MultiQC table
    columns_to_drop = [
        'gc_content', 'ambiguous_detail', 'n_regions', 
        'stop_positions', 'homopolymer_issues', 'frameshift_evidence'
    ]
    
    # Create copy and drop unwanted columns
    df_mqc = df.copy()
    for col in columns_to_drop:
        if col in df_mqc.columns:
            df_mqc = df_mqc.drop(columns=[col])
    
    # MultiQC table comment header
    table_comment = dedent("""
        # plot_type: 'table'
        # section_name: 'CSFV FASTA QC'
        # section_href: 'https://github.com/CFIA-NCFAD/nf-ionampliseq'
        # description: 'Quality metrics for CSFV consensus FASTA sequences.'
        # pconfig:
        #     namespace: 'CSFV'
        # headers:
        #     sequence_id:
        #         title: 'Sequence ID'
        #         description: 'Sequence identifier'
        #         scale: False
        #         format: '{}'
        #     quality_comment:
        #         title: 'QC Comment'
        #         description: 'Quality assessment summary'
        #         scale: False
        #         format: '{}'
        #     frameshift_evidence:
        #         title: 'Frameshift Evidence'
        #         description: 'Evidence for potential frameshifts'
        #         scale: False
        #         format: '{}'
        #     coverage_percentage:
        #         title: 'Genome Coverage (%)'
        #         description: 'Percentage of genome covered (excluding Ns)'
        #         format: '{:.1f}'
        #         min: 0
        #         max: 100
        #     length:
        #         title: 'Sequence Length'
        #         description: 'Total sequence length in nucleotides'
        #         format: '{:.0f}'
        #     total_ns:
        #         title: 'Total Ns'
        #         description: 'Total number of N characters in sequence'
        #         format: '{:.0f}'
        #     end_5prime_ns:
        #         title: 'End 5-prime Ns'
        #         description: 'Number of N characters at 5-prime end'
        #         format: '{:.0f}'
        #     end_3prime_ns:
        #         title: 'End 3-prime Ns'
        #         description: 'Number of N characters at 3-prime end'
        #         format: '{:.0f}'
        #     ambiguous_bases:
        #         title: 'Ambiguous Bases'
        #         description: 'Total number of ambiguous bases (not ATGCN)'
        #         format: '{:.0f}'
        #     polyprotein_start:
        #         title: 'Polyprotein Start'
        #         description: 'Start position of polyprotein region'
        #         format: '{:.0f}'
        #     polyprotein_end:
        #         title: 'Polyprotein End'
        #         description: 'End position of polyprotein region'
        #         format: '{:.0f}'
        #     potential_frameshift:
        #         title: 'Potential Frameshift'
        #         description: 'Whether potential frameshift was detected'
        #         scale: False
        #         format: '{}'
        #     stop_codons:
        #         title: 'Stop Codons'
        #         description: 'Number of premature stop codons in polyprotein'
        #         format: '{:.0f}'
        #     best_alignment_edit_distance:
        #         title: 'Best Alignment Edit Distance'
        #         description: 'Edit distance from reference polyprotein alignment using Edlib'
        #         format: '{:.0f}'
        #     best_alignment_score:
        #         title: 'Best Alignment Score'
        #         description: 'Alignment score from reference polyprotein alignment using Edlib'
        #         format: '{:.0f}'
        """)
    
    # Write MultiQC table with comment header
    with open(output_file, 'w') as fh:
        fh.write(table_comment)
        df_mqc.to_csv(fh, sep="\t", index=False)


def version_callback(value: bool):
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()

@app.command()
def analyze(
    input_dir: Path = typer.Argument(..., help="Directory containing FASTA files"),
    output_file: Path = typer.Option(
        "csfv_quality_summary.csv", "--output-csv", "-o", help="Output CSV filename"
    ),
    output_mqc: Path = typer.Option(
        "csfv_quality_summary_mqc.txt", "--output-mqc", "-m", help="Output MultiQC table filename"
    ),
    pattern: str = typer.Option(
        "*.fasta", help="File pattern to match (e.g., '*.fasta', '*.fa')"
    ),
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        count=True,
        help="Increase verbosity (-v for DEBUG, -vv for more verbose)",
    ),
    version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
) -> None:
    """Analyze CSFV FASTA consensus files for quality metrics."""

    set_log_level(verbose)

    # Validate input directory
    if not input_dir.exists():
        console.print(f"[red]Error: Directory {input_dir} does not exist[/red]")
        raise typer.Exit(1)

    if not input_dir.is_dir():
        console.print(f"[red]Error: {input_dir} is not a directory[/red]")
        raise typer.Exit(1)

    # Find FASTA files
    fasta_files = list(input_dir.glob(pattern))

    if not fasta_files:
        console.print(
            f"[red]No FASTA files found in {input_dir} matching pattern '{pattern}'[/red]"
        )
        console.print("Try patterns like: '*.fasta', '*.fa', '*.fas'")
        raise typer.Exit(1)

    console.print(f"[green]Found {len(fasta_files)} FASTA files[/green]")

    # Analyze files with progress bar
    all_results = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Analyzing files...", total=len(fasta_files))

        for fasta_file in fasta_files:
            logger.debug(f"Analyzing: {fasta_file}")
            progress.update(task, description=f"Analyzing {fasta_file.name}")

            try:
                results = analyze_fasta_file(fasta_file)
                all_results.extend(results)
                progress.advance(task)
            except Exception as e:
                logger.error(f"Failed to analyze {fasta_file}: {e}")
                continue

    if not all_results:
        console.print("[red]No sequences were successfully analyzed[/red]")
        raise typer.Exit(1)

    # Calculate average GC content for comparison
    avg_gc = statistics.mean(
        [r["gc_content"] for r in all_results if r["gc_content"] > 0]
    )

    # Convert dictionary fields to strings for CSV output
    for result in all_results:
        if isinstance(result.get("ambiguous_detail"), dict):
            result["ambiguous_detail"] = str(result["ambiguous_detail"])

    # Create DataFrame
    df = pd.DataFrame(all_results)

    # Sort alphanumerically by sequence_id - simple and clean
    df = df.sort_values("sequence_id")

    # Display summary statistics
    summary_table = Table(title="Analysis Summary")
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="green")

    summary_table.add_row("Total sequences", str(len(df)))
    summary_table.add_row("Average length", f"{df['length'].mean():.0f} bp")
    summary_table.add_row(
        "Average coverage", f"{df['coverage_percentage'].mean():.1f}%"
    )
    summary_table.add_row("Average GC content", f"{avg_gc:.1f}%")
    summary_table.add_row(
        "Sequences >95% coverage", str((df["coverage_percentage"] > 95).sum())
    )
    summary_table.add_row(
        "Sequences >90% coverage", str((df["coverage_percentage"] > 90).sum())
    )
    summary_table.add_row(
        "Sequences with potential frameshifts", str(df["potential_frameshift"].sum())
    )

    # Define column order for CSV output
    df_csv = df[[x for x, _ in COLUMNS]]
    df_csv = df_csv.rename(columns=dict(COLUMNS))  # type: ignore

    # Save CSV results
    df_csv.to_csv(output_file, index=False)
    console.print(f"[green]Results saved to: {output_file}[/green]")

    # Generate MultiQC table if requested
    if output_mqc:
        write_multiqc_table(df, output_mqc)
        console.print(f"[green]MultiQC table saved to: {output_mqc}[/green]")

    console.print(summary_table)


if __name__ == "__main__":
    app()

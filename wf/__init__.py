from typing import List, Optional

from latch.resources.workflow import workflow
from latch.types.directory import LatchDir
from latch.types.file import LatchFile

from wf.Parameters import metadata
from wf.task import *


@workflow(metadata)
def crispresso2_test(
    sample_name: str,
    run_name: str,
    output_folder: LatchDir,
    fastq_r1: LatchFile,
    amplicon_seq: List[str],
    fastq_r2: Optional[LatchFile] = None,
    guide_seq: Optional[List[str]] = None,
    amplicon_name: Optional[List[str]] = None,
    expected_hdr_amplicon_seq: Optional[str] = None,
    coding_seq: Optional[str] = None,
    guide_name: Optional[List[str]] = None,
    flexiguide_seq: Optional[str] = None,
    flexiguide_homology: Optional[int] = 80,
    flexiguide_name: Optional[str] = None,
    fastp_command: Optional[str] = "fastp",
    fastp_options_string: Optional[str] = None,
    quantification_window_coordinates: Optional[str] = None,
    amplicon_min_alignment_score: Optional[List[int]] = None,
    prime_editing_pegRNA_spacer_seq: Optional[str] = None,
    prime_editing_pegRNA_extension_seq: Optional[str] = None,
    prime_editing_pegRNA_scaffold_seq: Optional[str] = None,
    prime_editing_pegRNA_scaffold_min_match_length: Optional[int] = None,
    prime_editing_nicking_guide_seq: Optional[str] = None,
    prime_editing_override_prime_edited_ref_seq: Optional[str] = None,
    annotate_wildtype_allele: Optional[str] = None,
    file_prefix: Optional[str] = None,
    name: Optional[str] = None,
    dsODN: Optional[str] = None,
    bam_input: Optional[str] = None,
    bam_chr_loc: Optional[str] = None,
    custom_default_min_aln_score: Optional[int] = None,
    min_frequency_alleles_around_cut_to_plot: float = 0.2,
    max_rows_alleles_around_cut_to_plot: int = 50,
    exclude_bp_from_left: int = 15,
    exclude_bp_from_right: int = 15,
    min_average_read_quality: int = 0,
    min_single_bp_quality: int = 0,
    min_bp_quality_or_N: int = 0,
    conversion_nuc_from: str = "C",
    conversion_nuc_to: str = "T",
    prime_editing_pegRNA_extension_quantification_window_size: int = 5,
    quantification_window_size: int = 1,
    quantification_window_center: int = -3,
    default_min_aln_score: int = 60,
    plot_window_size: int = 20,
    discard_guide_positions_overhanging_amplicon_edge: bool = False,
    split_interleaved_input: bool = False,
    trim_sequences: bool = False,
    ignore_substitutions: bool = False,
    ignore_insertions: bool = False,
    ignore_deletions: bool = False,
    discard_indel_reads: bool = False,
    expand_ambiguous_alignments: bool = False,
    base_editor_output: bool = False,
    plot_histogram_outliers: bool = False,
    expand_allele_plots_by_quantification: bool = False,
    allele_plot_pcts_only_for_assigned_reference: bool = False,
    write_detailed_allele_table: bool = False,
    fastq_output: bool = False,
    keep_intermediate: bool = False,
    dump: bool = False,
    crispresso1_mode: bool = False,
    suppress_report: bool = False,
    suppress_plots: bool = False,
    place_report_in_output_folder: bool = False,
    debug: bool = False,
    no_rerun: bool = False,
    auto: bool = False,
    needleman_wunsch_aln_matrix_loc: str = "EDNAFULL",
    needleman_wunsch_gap_open: int = -20,
    needleman_wunsch_gap_extend: int = -2,
    needleman_wunsch_gap_incentive: int = 1,
    assign_ambiguous_alignments_to_first_reference: bool = False,
    suppress_amplicon_name_truncation: bool = False,
    verbosity: int = 3,
    zip_output: bool = False,
    disable_guardrails: bool = False,
) -> LatchDir:
    """Analysis of deep sequencing data for rapid and intuitive interpretation of genome editing experiments

    # CRISPResso2

    CRISPResso2 is a software pipeline designed to enable rapid and intuitive interpretation of genome editing experiments. A limited web implementation is available at: https://crispresso2.pinellolab.org/.

    Briefly, CRISPResso2:
    - aligns sequencing reads to a reference sequence
    - quantifies insertions, mutations and deletions to determine whether a read is modified or unmodified by genome editing
    - summarizes editing results in intuitive plots and datasets

    ## How to run CRISPResso2 on Latch ☕

    Check out our docs for detailed information on how to run CRISPResso2, including a how-to, and specifications on inputs and outputs: https://www.latch.wiki/crispresso2

    ## What can I do with CRISPResso2?

    CRISPResso2 can be used to analyze genome editing outcomes using cleaving nucleases (e.g. Cas9 or Cpf1) or noncleaving nucleases (e.g. base editors). The following operations can be automatically performed:
    - filtering of low-quality reads
    - adapter trimming
    - alignment of reads to one or multiple reference sequences (in the case of multiple alleles)
    - quantification of HDR and NHEJ outcomes (if the HDR sequence is provided)
    - quantification frameshift/inframe mutations and identification affected splice sites (if an exon sequence is provided)
    - visualization of the indel distribution and position (for cleaving nucleases)
    - visualization of distribution and position of substitutions (for base editors)
    - visualization of alleles and their frequencies"""

    run_dir = Initialize(outdir=output_folder, run_name=run_name)
    results = crispresso2(
        output_folder=run_dir,
        fastq_r1=fastq_r1,
        fastq_r2=fastq_r2,
        amplicon_seq=amplicon_seq,
        amplicon_name=amplicon_name,
        guide_seq=guide_seq,
        expected_hdr_amplicon_seq=expected_hdr_amplicon_seq,
        coding_seq=coding_seq,
        guide_name=guide_name,
        flexiguide_seq=flexiguide_seq,
        flexiguide_homology=flexiguide_homology,
        flexiguide_name=flexiguide_name,
        discard_guide_positions_overhanging_amplicon_edge=discard_guide_positions_overhanging_amplicon_edge,
        split_interleaved_input=split_interleaved_input,
        min_average_read_quality=min_average_read_quality,
        min_single_bp_quality=min_single_bp_quality,
        min_bp_quality_or_N=min_bp_quality_or_N,
        trim_sequences=trim_sequences,
        fastp_command=fastp_command,
        fastp_options_string=fastp_options_string,
        quantification_window_size=quantification_window_size,
        quantification_window_center=quantification_window_center,
        quantification_window_coordinates=quantification_window_coordinates,
        exclude_bp_from_left=exclude_bp_from_left,
        exclude_bp_from_right=exclude_bp_from_right,
        ignore_substitutions=ignore_substitutions,
        ignore_insertions=ignore_insertions,
        ignore_deletions=ignore_deletions,
        discard_indel_reads=discard_indel_reads,
        amplicon_min_alignment_score=amplicon_min_alignment_score,
        default_min_aln_score=default_min_aln_score,
        custom_default_min_aln_score=custom_default_min_aln_score,
        expand_ambiguous_alignments=expand_ambiguous_alignments,
        needleman_wunsch_gap_open=needleman_wunsch_gap_open,
        needleman_wunsch_gap_extend=needleman_wunsch_gap_extend,
        needleman_wunsch_gap_incentive=needleman_wunsch_gap_incentive,
        needleman_wunsch_aln_matrix_loc=needleman_wunsch_aln_matrix_loc,
        base_editor_output=base_editor_output,
        conversion_nuc_from=conversion_nuc_from,
        conversion_nuc_to=conversion_nuc_to,
        prime_editing_pegRNA_spacer_seq=prime_editing_pegRNA_spacer_seq,
        prime_editing_pegRNA_extension_seq=prime_editing_pegRNA_extension_seq,
        prime_editing_pegRNA_extension_quantification_window_size=prime_editing_pegRNA_extension_quantification_window_size,
        prime_editing_pegRNA_scaffold_seq=prime_editing_pegRNA_scaffold_seq,
        prime_editing_pegRNA_scaffold_min_match_length=prime_editing_pegRNA_scaffold_min_match_length,
        prime_editing_nicking_guide_seq=prime_editing_nicking_guide_seq,
        prime_editing_override_prime_edited_ref_seq=prime_editing_override_prime_edited_ref_seq,
        plot_histogram_outliers=plot_histogram_outliers,
        plot_window_size=plot_window_size,
        min_frequency_alleles_around_cut_to_plot=min_frequency_alleles_around_cut_to_plot,
        max_rows_alleles_around_cut_to_plot=max_rows_alleles_around_cut_to_plot,
        expand_allele_plots_by_quantification=expand_allele_plots_by_quantification,
        allele_plot_pcts_only_for_assigned_reference=allele_plot_pcts_only_for_assigned_reference,
        annotate_wildtype_allele=annotate_wildtype_allele,
        file_prefix=file_prefix,
        name=sample_name,
        write_detailed_allele_table=write_detailed_allele_table,
        fastq_output=fastq_output,
        keep_intermediate=keep_intermediate,
        dump=dump,
        crispresso1_mode=crispresso1_mode,
        suppress_report=suppress_report,
        suppress_plots=suppress_plots,
        place_report_in_output_folder=place_report_in_output_folder,
        auto=auto,
        dsODN=dsODN,
        debug=debug,
        no_rerun=no_rerun,
        bam_input=bam_input,
        bam_chr_loc=bam_chr_loc,
        assign_ambiguous_alignments_to_first_reference=assign_ambiguous_alignments_to_first_reference,
        suppress_amplicon_name_truncation=suppress_amplicon_name_truncation,
        verbosity=verbosity,
        zip_output=zip_output,
        disable_guardrails=disable_guardrails,
    )
    print(results)

    table_id = Update_Registry_Tables(
        outdir=results, run_name=run_name, sample=sample_name
    )
    return results

from latch.types import (
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    Params,
    Section,
    Spoiler,
)
from latch.types.metadata import (
    Multiselect,
    MultiselectOption,
)

flow = [
    Section(
        "Inputs/Outputs",
        Params(
            "fastq_r1",
            "fastq_r2",
            "sample_name",
            "run_name",
            "output_folder",
        ),
        Spoiler(
            "Optional",
            Params(
                "split_interleaved_input",
            ),
        ),
    ),
    Section(
        "Amplicon Sequence Parameters",
        Params("amplicon_seq"),
        Spoiler(
            "Optional",
            Params(
                "amplicon_name",
                "expected_hdr_amplicon_seq",
                "default_min_aln_score",
                "suppress_amplicon_name_truncation",
            ),
        ),
    ),
    Section(
        "sgRNA settings",
        Params("guide_seq"),
        Spoiler(
            "Optional",
            Params(
                "",
                "quantification_window_size",
                "quantification_window_center",
                "discard_guide_positions_overhanging_amplicon_edge",
                "ignore_substitutions",
                "ignore_insertions",
                "ignore_deletions",
                "discard_indel_reads",
                "plot_window_size",
            ),
        ),
    ),
    Spoiler(
        "Optional Crispresso2 Arguments",
        Spoiler(
            "FastP Options",
            Params("trim_sequences", "fastp_command", "fastp_options_string"),
        ),
        Spoiler(
            "Prime Editing Parameters",
            Params(
                "prime_editing_pegRNA_extension_quantification_window_size",
                "prime_editing_pegRNA_spacer_seq" "prime_editing_pegRNA_extension_seq",
                "prime_editing_pegRNA_scaffold_seq",
                "prime_editing_pegRNA_scaffold_min_match_length",
                "prime_editing_nicking_guide_seq",
                "prime_editing_override_prime_edited_ref_seq",
            ),
        ),
        Spoiler(
            "Base Editing Parameters",
            Params(
                "base_editor_output",
                "conversion_nuc_from",
                "conversion_nuc_to",
                "dsODN",
            ),
        ),
        Spoiler(
            "Flexiguide Parameters",
            Params("flexiguide_seq", "flexiguide_homology", "flexiguide_name"),
        ),
        Spoiler(
            "Allele Plot Parameters",
            Params(
                "annotate_wildtype_allele",
                "min_frequency_alleles_around_cut_to_plot",
                "max_rows_alleles_around_cut_to_plot",
                "expand_allele_plots_by_quantification",
                "allele_plot_pcts_only_for_assigned_reference",
                "write_detailed_allele_table",
            ),
        ),
        Spoiler(
            "Alignment",
            Params(
                "expand_ambiguous_alignments",
                "needleman_wunsch_aln_matrix_loc",
                "needleman_wunsch_gap_open",
                "needleman_wunsch_gap_extend",
                "needleman_wunsch_gap_incentive",
                "assign_ambiguous_alignments_to_first_reference",
            ),
        ),
        Spoiler(
            "Quality Filtering",
            Params(
                "min_average_read_quality",
                "min_single_bp_quality",
                "min_bp_quality_or_N",
            ),
        ),
        Spoiler(
            "base pair exclusion from amplicon sequence for quantification of mutations",
            Params("exclude_bp_from_left", "exclude_bp_from_right"),
        ),
        Spoiler(
            "Output settings",
            Params(
                "plot_histogram_outliers",
                "fastq_output",
                "crispresso1_mode",
                "place_report_in_output_folder",
                "file_prefix",
                "supress_plots",
                "supress_report",
                "zip_output",
            ),
        ),
        Spoiler("Bam instead of FastQ", Params("bam_input", "bam_chr_loc")),
        Spoiler(
            "Debug Parameters",
            Params(
                "keep_intermediate",
                "dump",
                "debug",
                "no_rerun",
                "auto",
                "verbosity",
                "disable_guardrails",
            ),
        ),
    ),
]


metadata = LatchMetadata(
    display_name="Crispresso2",
    author=LatchAuthor(
        name="Harihara",
    ),
    parameters={
        "fastq_r1": LatchParameter(
            display_name="Fastq R1",
            description="The first FastQ file, if you use only a fastq file here CRISPResso will assume it contains single end reads and run it as such.",
        ),
        "split_interleaved_input": LatchParameter(
            description="Marks the file in Read 1 as containing interleaved reads. CRISPResso will split the paired end reads into two files before running.",
            display_name="Read 1 is Interleaved",
        ),
        "suppress_amplicon_name_truncation": LatchParameter(
            display_name="Suppress amplicon name truncation",
            description="If set, amplicon names will not be truncated when creating output filename prefixes. If not set, amplicon names longer than 21 characters will be truncated when creating filename prefixes.",
        ),
        "fastq_r2": LatchParameter(
            display_name="Fastq R2",
            description="The second fastq file for paired end reads.",
        ),
        "amplicon_seq": LatchParameter(
            description="The amplicon sequence(s) used for the experiment.",
            display_name="Amplicon Sequence",
        ),
        "amplicon_name": LatchParameter(
            description="A name for each reference amplicon can be given. If multiple amplicons are given, multiple names can be specified here and their order must correspond to the amplicons given.",
            display_name="Amplicon Name",
        ),
        "guide_seq": LatchParameter(
            description="sgRNAs should be input as the guide RNA sequence (usually 20 nt) immediately adjacent to but not including the PAM sequence (5' (left) of NGG for SpCas9). If the sgRNA is not provided, quantification may include modifications far from the predicted editing site and may result in overestimation of editing rates.",
            display_name="Guide Sequence",
        ),
        "guide_name": LatchParameter(
            description="A name for each sgRNA can be given. If multiple sgRNA are given, multiple names can be specified here and their order must correspond to the sgRNA given.",
            display_name="Guide Sequence Name",
        ),
        "output_folder": LatchParameter(
            description="A file path pointing to where your results will live.",
            display_name="Output Location",
        ),
        "sample_name": LatchParameter(display_name="Sample Name"),
        "run_name": LatchParameter(display_name="Run Name"),
        "default_min_aln_score": LatchParameter(
            description="Sequences must have at least this homology percentage score with the amplicon to be aligned.",
            display_name="Minimum Homology % For Alignment to an Amplicon",
            appearance_type=Multiselect([50, 60, 70, 80, 90], allow_custom=True),
        ),
        "discard_guide_positions_overhanging_amplicon_edge": LatchParameter(
            description="If set, for guides that align to multiple positions, guide positions will be discarded if plotting around those regions would included bp that extend beyond the end of the amplicon.",
            display_name="Discard Guide Positions That Extend Beyond End of Amplicon",
        ),
        "quantification_window_center": LatchParameter(
            description="Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17.",
            display_name="Quantification Window Center",
            appearance_type=Multiselect(
                [
                    MultiselectOption("Cas9", -3),
                    MultiselectOption("Cpfl1", 1),
                    MultiselectOption("Base Editors", -10),
                ],
                allow_custom=True,
            ),
        ),
        "assign_ambiguous_alignments_to_first_reference": LatchParameter(
            display_name="Assign ambiguous alignments to first reference",
            description="If more than one reference amplicon is given, ambiguous reads that align with the same score to multiple amplicons will be assigned to the first amplicon. Default behavior is to exclude ambiguous alignments.",
        ),
        "quantification_window_size": LatchParameter(
            description='Defines the size (in bp) of the quantification window extending from the position specified by the "--cleavage_offset" or "--quantification_window_center" parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp.',
            display_name="Quantification Window Size",
            appearance_type=Multiselect(
                [
                    MultiselectOption("Cas9", -3),
                    MultiselectOption("Cpfl1", 1),
                    MultiselectOption("Base Editors", -10),
                    MultiselectOption("Disable", 0),
                ],
                allow_custom=True,
            ),
        ),
        "plot_window_size": LatchParameter(
            description="Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are plotted.",
            display_name="Plot Window Size",
            appearance_type=Multiselect([5, 10, 20, 30], allow_custom=True),
        ),
        "ignore_substitutions": LatchParameter(
            description="Ignore substitutions events for the quantification and visualization.",
            display_name="Ignore Substitutions",
        ),
        "ignore_insertions": LatchParameter(
            description="Ignore insertions events for the quantification and visualization.",
            display_name="Ignore Insertions",
        ),
        "ignore_deletions": LatchParameter(
            description="Ignore deletions events for the quantification and visualization.",
            display_name="Ignore Deletions",
        ),
        "discard_indel_reads": LatchParameter(
            description="Discard reads with indels in the quantification window from analysis.",
            display_name="Ignore Indel Reads",
        ),
        "prime_editing_pegRNA_spacer_seq": LatchParameter(
            description="pegRNA spacer sgRNA sequence used in prime editing. The spacer should not include the PAM sequence. The sequence should be given in the RNA 5'->3' order, so for Cas9, the PAM would be on the right side of the given sequence.",
            display_name="Spacer Sequence",
        ),
        "prime_editing_pegRNA_extension_seq": LatchParameter(
            description="Extension sequence used in prime editing. The sequence should be given in the RNA 5'->3' order, such that the sequence starts with the RT template including the edit, followed by the Primer-binding site (PBS).",
            display_name="Extension Sequence",
        ),
        "prime_editing_pegRNA_extension_quantification_window_size": LatchParameter(
            description="Quantification window size (in bp) at flap site for measuring modifications anchored at the right side of the extension sequence. Similar to the --quantification_window parameter, the total length of the quantification window will be 2x this parameter. Default: 5bp (10bp total window size)",
            display_name="pegRNA Extension Quantification Window Size",
            appearance_type=Multiselect([1, 5, 10], allow_custom=True),
        ),
        "prime_editing_nicking_guide_seq": LatchParameter(
            description="Nicking sgRNA sequence used in prime editing. The sgRNA should not include the PAM sequence. The sequence should be given in the RNA 5'->3' order, so for Cas9, the PAM would be on the right side of the sequence.",
            display_name="Nicking sgRNA",
        ),
        "prime_editing_pegRNA_scaffold_seq": LatchParameter(
            description="If given, reads containing any of this scaffold sequence before extension sequence (provided by --prime_editing_extension_seq) will be classified as 'Scaffold-incorporated'. The sequence should be given in the 5'->3' order such that the RT template directly follows this sequence. A common value ends with 'GGCACCGAGUCGGUGC'.",
            display_name="Scaffold Sequence",
        ),
        "prime_editing_pegRNA_scaffold_min_match_length": LatchParameter(
            description="Minimum number of bases matching scaffold sequence for the read to be counted as 'Scaffold-incorporated'. If the scaffold sequence matches the reference sequence at the incorporation site, the minimum number of bases to match will be minimally increased (beyond this parameter) to disambiguate between prime-edited and scaffold-incorporated sequences.",
            display_name="Scaffold Sequence",
        ),
        "prime_editing_override_prime_edited_ref_seq": LatchParameter(
            description="If given, this sequence will be used as the prime-edited reference sequence. This may be useful if the prime-edited reference sequence has large indels or the algorithm cannot otherwise infer the correct reference sequence.",
            display_name="Specify Prime Edited Reference Sequence",
        ),
        "base_editor_output": LatchParameter(
            description="Outputs plots and tables to aid in analysis of base editor studies. If base editor output is selected, plots showing the frequency of substitutions in the quantification window are generated. The target and result bases can also be set to measure the rate of on-target conversion at bases in the quantification window.",
            display_name="Base Editing Output",
        ),
        "conversion_nuc_from": LatchParameter(
            description="For base editor plots, this is the nucleotide targeted by the base editor.",
            display_name="Base Editor Target Base",
        ),
        "conversion_nuc_to": LatchParameter(
            description="For base editor plots, this is the nucleotide produced by the base editor.",
            display_name="Base Editor Output Base",
        ),
        "expected_hdr_amplicon_seq": LatchParameter(
            description="Amplicon sequence expected after HDR. The expected HDR amplicon sequence can be provided to quantify the number of reads showing a successful HDR repair.",
            display_name="Include HDR Sequence",
        ),
        "coding_seq": LatchParameter(
            description="Subsequence/s of the amplicon sequence covering one or more coding sequences for frameshift analysis. Sequences of exons within the amplicon sequence can be provided to enable frameshift analysis and splice site analysis by CRISPResso2. If more than one (for example, split by intron/s), please separate by commas. Users should provide the subsequences of the reference amplicon sequence that correspond to coding sequences (not the whole exon sequence(s)!).",
            display_name="Include Exon Coding Sequence",
        ),
        "dsODN": LatchParameter(
            description="dsODN sequence -- Reads containing the dsODN are labeled and quantified.",
            display_name="Include dsODN Sequence",
        ),
        "min_average_read_quality": LatchParameter(
            description="Minimum average quality score (phred33) to keep a read.",
            display_name="Minimum Average Read Quality",
            appearance_type=Multiselect(
                [
                    MultiselectOption("No Filter", 0),
                    MultiselectOption(">10", 10),
                    MultiselectOption(">20", 20),
                    MultiselectOption(">30", 30),
                    MultiselectOption(">35", 35),
                ]
            ),
        ),
        "min_single_bp_quality": LatchParameter(
            description="Minimum single bp score (phred33) to keep a read (default: 0)",
            display_name="Minimum Single Base Pair Quality",
            appearance_type=Multiselect(
                [
                    MultiselectOption("No Filter", 0),
                    MultiselectOption(">10", 10),
                    MultiselectOption(">20", 20),
                    MultiselectOption(">30", 30),
                    MultiselectOption(">35", 35),
                ]
            ),
        ),
        "min_bp_quality_or_N": LatchParameter(
            description='Bases with a quality score (phred33) less than this value will be set to "N".',
            display_name="Replace Bases With N That Have a Quality Lower Than",
            appearance_type=Multiselect(
                [
                    MultiselectOption("No Filter", 0),
                    MultiselectOption(">10", 10),
                    MultiselectOption(">20", 20),
                    MultiselectOption(">30", 30),
                    MultiselectOption(">35", 35),
                ]
            ),
        ),
        "exclude_bp_from_left": LatchParameter(
            description="Exclude bp from the left side of the amplicon sequence for the quantification of the indels.",
            display_name="Base Pairs Excluded from the Left Side",
            appearance_type=Multiselect(
                [
                    MultiselectOption("Disabled", 0),
                    MultiselectOption("5", 5),
                    MultiselectOption("15", 15),
                    MultiselectOption("20", 20),
                    MultiselectOption("40", 40),
                    MultiselectOption("50", 50),
                ],
                allow_custom=True,
            ),
        ),
        "exclude_bp_from_right": LatchParameter(
            description="Exclude bp from the right side of the amplicon sequence for the quantification of the indels.",
            display_name="Base Pairs Excluded from the Right Side",
            appearance_type=Multiselect(
                [
                    MultiselectOption("Disabled", 0),
                    MultiselectOption("5", 5),
                    MultiselectOption("15", 15),
                    MultiselectOption("20", 20),
                    MultiselectOption("40", 40),
                    MultiselectOption("50", 50),
                ],
                allow_custom=True,
            ),
        ),
        "trim_sequences": LatchParameter(
            description="Enable the trimming of Illumina adapters with FastP.",
            display_name="Enable Trimming Adaptor",
        ),
        "fastp_options_string": LatchParameter(
            description="Override options for fastp, e.g. --length_required 70 --umi",
            display_name="FastP Options String",
        ),
        "expand_ambiguous_alignments": LatchParameter(
            description="If more than one reference amplicon is given, reads that align to multiple reference amplicons will count equally toward each amplicon. Default behavior is to exclude ambiguous alignments.",
            display_name="Expand Ambiguous Alignments",
        ),
        "needleman_wunsch_gap_open": LatchParameter(
            description="Gap open option for Needleman-Wunsch alignment.",
            display_name="Needleman Wunsch Gap Open",
        ),
        "needleman_wunsch_gap_incentive": LatchParameter(
            description="Gap incentive value for inserting indels at cut sites.",
            display_name="Needleman Wunsch Gap Incentive",
        ),
        "needleman_wunsch_gap_extend": LatchParameter(
            description="Gap extend option for Needleman-Wunsch alignment.",
            display_name="Needleman Wunsch Gap Extend",
        ),
        "plot_histogram_outliers": LatchParameter(
            description="If set, all values will be shown on histograms. By default (if unset), histogram ranges are limited to plotting data within the 99 percentile.",
            display_name="Plot Histogram Outliers",
        ),
        "max_rows_alleles_around_cut_to_plot": LatchParameter(
            description="Maximum number of rows to report in the alleles table plot.",
            display_name="Maximum Rows Reported in the Alleles Table",
            appearance_type=Multiselect([5, 10, 20, 30, 40, 50], allow_custom=True),
        ),
        "min_frequency_alleles_around_cut_to_plot": LatchParameter(
            description="Minimum % reads required to report an allele in the alleles table plot. This parameter only affects plotting. All alleles will be reported in data files.",
            display_name="Minimum % Reads Required To Report an Allele in Table Plot",
            appearance_type=Multiselect([0.2, 0.4, 0.6, 0.8], allow_custom=True),
        ),
        "annotate_wildtype_allele": LatchParameter(
            description="Wildtype alleles in the allele table plots will be marked with this string (e.g. **).",
            display_name="Annotate Wildtype Allele",
        ),
        "allele_plot_pcts_only_for_assigned_reference": LatchParameter(
            description="If set, in the allele plots, the percentages will show the percentage as a percent of reads aligned to the assigned reference. Default behavior is to show percentage as a percent of all reads.",
            display_name="Show Percentage As Reads Aligned to Assigned Reference",
        ),
        "write_detailed_allele_table": LatchParameter(
            description="If set, a detailed allele table will be written including alignment scores for each read sequence.",
            display_name="Include Alignment Scores for Each Read Sequence in Allele Table",
        ),
        "expand_allele_plots_by_quantification": LatchParameter(
            description="If set, alleles with different modifications in the quantification window (but not necessarily in the plotting window (e.g. for another sgRNA)) are plotted on separate lines, even though they may have the same apparent sequence. To force the allele plot and the allele table to be the same, set this parameter. If unset, all alleles with the same sequence will be collapsed into one row.",
            display_name="Force Allele Plot and Allele Table To Be the Same",
        ),
        "fastq_output": LatchParameter(
            description="If set, a fastq file with annotations for each read will be produced.",
            display_name="Output FastQ File for Each Read",
        ),
        "crispresso1_mode": LatchParameter(
            description="Output as in CRISPResso1. In particular, if this flag is set, the old output files 'Mapping_statistics.txt', and 'Quantification_of_editing_frequency.txt' are created, and the new files 'nucleotide_frequency_table.txt' and 'substitution_frequency_table.txt' and figure 2a and 2b are suppressed, and the files 'selected_nucleotide_percentage_table.txt' are not produced when the flag --base_editor_output is set.",
            display_name="use CRISPResso1 Output Mode",
        ),
        "place_report_in_output_folder": LatchParameter(
            description="If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the report output.",
            display_name="Place Report in Same Folder as Output Data",
        ),
        "file_prefix": LatchParameter(
            description="File prefix for output plots and table.",
            display_name="Set File Prefix For Plots and Tables",
        ),
        "suppress_report": LatchParameter(
            description="Suppress output report, plots output as .pdf only (not .png).",
            display_name="Plot Output as PDF Only",
        ),
        "suppress_plots": LatchParameter(
            description="Suppress output plots.", display_name="Suppress Output Plots"
        ),
        "bam_input": LatchParameter(
            description="Aligned reads for processing in bam format. This parameter can be given instead of fastq_r1 to specify that reads are to be taken from this bam file. An output bam is produced that contains an additional field with CRISPResso2 information.",
            display_name="Read 1 (BAM)",
        ),
        "bam_chr_loc": LatchParameter(
            description='Chromosome location in bam for reads to process. For example: "chr1:50-100" or "chrX".',
            display_name="BAM Chromosome Location",
        ),
        "flexiguide_seq": LatchParameter(
            display_name="Flexiguide Seq",
            description="sgRNA sequence (flexible) (can be comma-separated list of multiple flexiguides). The flexiguide sequence will be aligned to the amplicon sequence(s), as long as the guide sequence has homology as set by --flexiguide_homology.",
        ),
        "flexiguide_homology": LatchParameter(
            display_name="Flexiguide Homology",
            description="flexiguides will yield guides in amplicons with at least this homology to the flexiguide sequence.",
            appearance_type=Multiselect([60, 70, 80, 90, 95], allow_custom=True),
        ),
        "flexiguide_name": LatchParameter(display_name="Flexiguide Name"),
        "verbosity": LatchParameter(
            display_name="verbosity",
            description="Help: Verbosity level of output to the console (1-4) 4 is the most verbose",
            appearance_type=Multiselect([1, 2, 3, 4], allow_custom=False),
        ),
        "zip_output": LatchParameter(
            display_name="Zip Output",
            description="If set, the output will be placed in a zip folder.",
        ),
        "disable_guardrails": LatchParameter(
            display_name="Disable Guardrails", description="Disable guardrail warnings"
        ),
    },
    flow=flow,
)

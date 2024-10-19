import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd
from latch.account import Account
from latch.registry.project import Project
from latch.registry.table import Table
from latch.resources.tasks import small_task
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.Parameters import Guide_RNA


@small_task
def Initialize(outdir: LatchDir, run_name: str) -> LatchDir:
    if os.path.isdir(f"{outdir.local_path}/{run_name}"):
        return LatchDir(f"{outdir.remote_directory}/{run_name}")
    else:
        os.mkdir(run_name)
        return LatchDir(run_name, os.path.join(outdir.remote_directory, run_name))


@small_task
def crispresso2(
    output_folder: LatchDir,
    run_name: str,
    fastq_r1: LatchFile,
    amplicon_seq: List[str],
    sgRNA_Seq: List[Guide_RNA],
    fastq_r2: Optional[LatchFile] = None,
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
    quantification_window_size: int = -3,
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
    if name is not None:
        name = name.replace(" ", "_")

    amplicon_seq = ",".join(amplicon_seq)
    amplicon_seq = amplicon_seq.upper()

    if amplicon_name is not None:
        amplicon_name = ",".join(amplicon_name)

    if guide_name is not None:
        guide_name = ",".join(guide_name)
    if amplicon_min_alignment_score is not None:
        amplicon_min_alignment_score = ",".join(
            [str(x) for x in amplicon_min_alignment_score]
        )
    if prime_editing_pegRNA_spacer_seq != None:
        prime_editing_pegRNA_spacer_seq = prime_editing_pegRNA_spacer_seq.upper()
    if prime_editing_pegRNA_extension_seq != None:
        prime_editing_pegRNA_extension_seq = prime_editing_pegRNA_extension_seq.upper()
    if prime_editing_pegRNA_scaffold_seq != None:
        prime_editing_pegRNA_scaffold_seq = prime_editing_pegRNA_scaffold_seq.upper()
    if prime_editing_nicking_guide_seq != None:
        prime_editing_nicking_guide_seq = prime_editing_nicking_guide_seq.upper()
    if prime_editing_override_prime_edited_ref_seq != None:
        prime_editing_override_prime_edited_ref_seq = (
            prime_editing_override_prime_edited_ref_seq.upper()
        )
    if flexiguide_seq != None:
        flexiguide_seq = flexiguide_seq.upper()

    header_args = locals()

    flag_params = [
        "base_editor_output",
        "discard_guide_positions_overhanging_amplicon_edge",
        "split_interleaved_input",
        "trim_sequences",
        "stringent_flash_merging",
        "ignore_substitutions",
        "ignore_insertions",
        "ignore_deletions",
        "discard_indel_reads",
        "expand_ambiguous_alignments",
        "base_editor_output",
        "plot_histogram_outliers",
        "expand_allele_plots_by_quantification",
        "allele_plot_pcts_only_for_assigned_reference",
        "write_detailed_allele_table",
        "fastq_output",
        "keep_intermediate",
        "dump",
        "crispresso1_mode",
        "suppress_report",
        "suppress_plots",
        "place_report_in_output_folder",
        "auto",
        "debug",
        "no_rerun",
        "assign_ambiguous_alignments_to_first_reference",
        "suppress_amplicon_name_truncation",
        "zip_output",
        "disable_guardrails",
    ]

    protected_params = {
        "output_folder": True,
        "fastq_r1": True,
        "fastq_r2": True,
        "amplicon_seq": True,
        "name": True,
        "run_name": True,
        "sgRNA_Seq": True,
    }

    optional_cmd_params = []
    for param, value in header_args.items():
        if param not in protected_params and value is not None:
            if param not in flag_params:
                if str(value) is not "" and value is not None and value != "<nil>":
                    optional_cmd_params.extend((f"--{param}", str(value)))
            else:
                if value:
                    optional_cmd_params.append(f"--{param}")

    fastq_r1_path = fastq_r1.local_path
    if fastq_r2 is not None:
        _ = fastq_r2.local_path

    if name is None or name == "":
        name = Path(fastq_r1_path).name.split(".")[0]

    if sgRNA_Seq is not None:
        guide_seqs = []
        for c in sgRNA_Seq:
            guide_seqs.append(c.guide_seq)
        guide_rna = ",".join(guide_seqs)
        guide_rna = guide_rna.upper()

    local_output_dir = Path(f"/root/outputs")
    local_output_dir.mkdir(parents=True, exist_ok=True)

    local_run_dir = Path(f"/root/outputs/{run_name}")
    local_run_dir.mkdir(parents=True, exist_ok=True)
    crisp_output_dir = Path(f"/root/outputs/{run_name}/CRISPResso_on_{name}")
    crisp_output_dir.mkdir(parents=True, exist_ok=True)
    # if not os.path.isdir(crisp_output_dir):
    #    os.mkdir(crisp_output_dir)

    crispresso_cmd = [
        "mamba",
        "run",
        "-n",
        "crispresso2_env",
        "CRISPResso",
        "--fastq_r1",
        fastq_r1.local_path,
        "--amplicon_seq",
        amplicon_seq,
        "--output_folder",
        str(crisp_output_dir),
    ]
    if sgRNA_Seq is not None:
        crispresso_cmd += ["--guide_seq", guide_rna]
    crispresso_cmd += optional_cmd_params

    crispresso_cmd.extend(("--name", name))
    if fastq_r2 is not None:
        crispresso_cmd += ["--fastq_r2", fastq_r2.local_path]

    print(f"Command: {crispresso_cmd}")
    cmd = " ".join(crispresso_cmd)
    results = subprocess.run(cmd, shell=True, check=True)

    print(f"Working directory {os.getcwd()}")
    if results is not None:
        print(f" stdout: {results.stdout}" + f" stderr: {results.stderr}")
    else:
        print("No process output generated")

    return LatchOutputDir(
        str(local_output_dir.resolve()),
        output_folder.remote_directory,
    )


@small_task
def update_registry_tables(outdir: LatchOutputDir, run_name: str, sample: str) -> str:
    print(outdir)
    sample_dir = f"{outdir.remote_path}/{run_name}/CRISPResso_on_{sample}/CRISPResso_on_{sample}/"
    for g in LatchDir(sample_dir).iterdir():
        fname = Path(g.remote_path).name
        if fname.startswith("Alleles_frequency_table_around_sgRNA_"):
            allele_freq_file = g
    mapping_statistics_file = LatchFile(
        f"{sample_dir}CRISPResso_mapping_statistics.txt"
    )
    editing_quantification = LatchFile(
        f"{sample_dir}CRISPResso_quantification_of_editing_frequency.txt"
    )

    df_editing_quantification = pd.read_csv(editing_quantification, sep="\t")
    df = pd.read_csv(mapping_statistics_file.local_path, sep="\t").T.to_dict()[0]
    d = {}
    num_reads = df["READS IN INPUTS"]

    d["num_reads"] = df["READS IN INPUTS"]

    d["num_reads_after_preprocessing"] = df["READS AFTER PREPROCESSING"]
    d["frac_reads_after_preprocessing"] = (
        d["num_reads_after_preprocessing"] / num_reads * 100.0
    )
    d["num_reads_aligned"] = df["READS ALIGNED"]
    d["frac_reads_aligned"] = d["num_reads_aligned"] / num_reads * 100.0
    d["num_discarded_reads"] = num_reads - d["num_reads_aligned"]
    d["frac_discarded_reads"] = d["num_discarded_reads"] / num_reads * 100.0
    d["Unmodified%"] = float(df_editing_quantification.iloc[0]["Unmodified%"])
    d["Modified%"] = float(df_editing_quantification.iloc[0]["Modified%"])
    d["Deletions"] = float(df_editing_quantification.iloc[0]["Deletions"])
    d["Deletion%"] = d["Deletions"] / d["num_reads"] * 100.0

    columns = [
        "num_reads",
        "num_reads_after_preprocessing",
        "frac_reads_after_preprocessing",
        "num_reads_aligned",
        "frac_reads_aligned",
        "num_discarded_reads",
        "frac_discarded_reads",
        "Unmodified%",
        "Modified%",
        "Deletions",
        "Deletion%",
    ]
    col_types = [float] * 11

    try:
        df_indels = pd.read_csv(allele_freq_file.local_path, sep="\t")
        df_indels_only = df_indels[
            ((df_indels["n_deleted"] > 0) | (df_indels["n_inserted"] > 0))
        ]

        top_5_alleles = df_indels.head(5)
        top_5_alleles[["Aligned_Sequence", "%Reads"]]
        alleles = top_5_alleles["Aligned_Sequence"].tolist()
        read_fracs = top_5_alleles["%Reads"].tolist()

        for i in range(len(alleles)):
            d[f"Allele_{str(i+1)}"] = alleles[i]
            d[f"Read_Frac_Allele_{str(i+1)}"] = read_fracs[i]
            columns.append(f"Allele_{str(i+1)}")
            col_types.append(str)
            columns.append(f"Read_Frac_Allele_{str(i+1)}")
            col_types.append(float)
    except Exception:
        print("Allele Frequency Tables missing. Probably a prime editing experiment? ")

    columns = ["sample"] + columns
    col_types = [str] + col_types
    d["sample"] = sample
    print(d)
    target_project_name = "Crispresso2_Runs"
    target_table_name = f"{run_name}_Results"
    account = Account.current()
    target_project = next(
        (
            project
            for project in account.list_registry_projects()
            if project.get_display_name() == target_project_name
        ),
        None,
    )

    if target_project is None:
        with account.update() as account_updater:
            account_updater.upsert_registry_project(target_project_name)
        target_project = next(
            project
            for project in account.list_registry_projects()
            if project.get_display_name() == target_project_name
        )
        print("Upserted project")

    target_table = next(
        (
            table
            for table in target_project.list_tables()
            if table.get_display_name() == target_table_name
        ),
        None,
    )

    if target_table == None:
        ### Upsert_Table
        with target_project.update() as project_updater:
            project_updater.upsert_table(target_table_name)
        target_table = next(
            table
            for table in target_project.list_tables()
            if table.get_display_name() == target_table_name
        )
        print("Upserted table")

        with target_table.update() as updater:
            for i in range(0, len(columns)):
                updater.upsert_column(columns[i], type=col_types[i])
        print("Upserted columns")
    print(d)

    with target_table.update() as updater:
        ctr = len(target_table.get_dataframe())
        updater.upsert_record(name=str(ctr), **d)

    return target_table.id

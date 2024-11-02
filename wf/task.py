import os
import subprocess
from dataclasses import dataclass
from os import listdir, mkdir
from os.path import isdir, isfile
from pathlib import Path
from shutil import copy, move
from typing import Optional

from dataclasses_json import dataclass_json
from latch import custom_task, small_task, workflow
from latch.account import Account
from latch.registry.project import Project
from latch.registry.table import Table
from latch.types import LatchDir, LatchFile, LatchOutputDir

import wf.Compress_Coverages


def Copy_Directory(remote_dir: LatchDir, local_dir: str):
    mkdir(local_dir)
    iter_dir = remote_dir.iterdir()
    for i in iter_dir:
        copy(i.local_path, local_dir)


@small_task
def Initialize(
    output_dir: LatchDir, sample: str, run_name: str
) -> (LatchDir, LatchDir, LatchDir, str, str):
    os.mkdir(run_name)
    os.mkdir(f"{run_name}/{sample}")
    os.mkdir(f"{run_name}/{sample}/Replicate_1")
    os.mkdir(f"{run_name}/{sample}/Replicate_2")
    outdir = LatchDir(sample, os.path.join(output_dir.remote_path, run_name))
    return (
        LatchDir(sample, os.path.join(outdir.remote_path, sample)),
        LatchDir(
            sample + "/Replicate_1",
            os.path.join(outdir.remote_path, sample, "Replicate_1"),
        ),
        LatchDir(
            sample + "/Replicate_2",
            os.path.join(outdir.remote_path, sample, "Replicate_2"),
        ),
        sample + ".Replicate_1",
        sample + ".Replicate_2",
    )


@custom_task(cpu=8, memory=48, storage_gib=300)
def Trim_Reads(
    fastq_1: Optional[LatchFile],
    fastq_2: Optional[LatchFile],
    sample: str,
    output_dir: LatchDir,
) -> (LatchFile, LatchFile, LatchDir, LatchDir):
    local_directory = sample + "/"
    if not os.path.isdir(local_directory):
        os.mkdir(local_directory)

    trim_galore_command = [
        "~/TrimGalore-0.6.10/trim_galore",
        "--cores",
        "8",
        "--paired",
        fastq_1.local_path,
        fastq_2.local_path,
        "--fastqc",
        "--output_dir",
        local_directory,
    ]
    subprocess.run(" ".join(trim_galore_command), shell=True, check=True)

    filepath = Path(fastq_1)
    trimmed_fastq_1 = filepath.name.replace(".fastq.gz", "_val_1.fq.gz")

    filepath = Path(fastq_2)
    trimmed_fastq_2 = filepath.name.replace(".fastq.gz", "_val_2.fq.gz")

    return (
        LatchFile(
            local_directory + "/" + trimmed_fastq_1,
            output_dir.remote_directory + "/trimmed_reads/" + trimmed_fastq_1,
        ),
        LatchFile(
            local_directory + "/" + trimmed_fastq_2,
            output_dir.remote_directory + "/trimmed_reads/" + trimmed_fastq_2,
        ),
        LatchDir(local_directory, output_dir.remote_directory + "/trimmed_reads/"),
        output_dir,
    )


@custom_task(cpu=24, memory=48, storage_gib=300)
def Align_Reads(
    trimmed_fastq_1: LatchFile,
    trimmed_fastq_2: LatchFile,
    genome: LatchDir,
    sample: str,
    output_dir: LatchDir,
) -> (LatchFile, LatchDir):
    idx_path = "BowTie2_Index/"
    Copy_Directory(genome, idx_path)

    files = os.listdir(idx_path)

    for f in files:
        if f.endswith(".1.bt2") and "rev" not in f:
            idx_prefix = idx_path + f.replace(".1.bt2", "")
            break
    print(idx_prefix)

    # local_bam_file = Path(trimmed_fastq_1.local_path).name.replace(
    #    "_val_1.fq.gz", ".bam"
    # )
    local_bam_file = str(sample) + ".bam"

    cmd_align = [
        "mamba",
        "run",
        "-n",
        "bowtie2_env",
        "bowtie2",
        "-p",
        "24",
        "-X2000",
        "--local",
        "--no-mixed",
        "--no-discordant",
        "-x",
        idx_prefix,
        "-1",
        trimmed_fastq_1.local_path,
        "-2",
        trimmed_fastq_2.local_path,
        "|",
        "samtools",
        "view",
        "-bS",
        "-",
        ">",
        local_bam_file,
    ]
    print(" ".join(cmd_align))
    subprocess.run(" ".join(cmd_align), shell=True, check=True)

    return LatchFile(
        Path(local_bam_file).resolve(),
        output_dir.remote_directory + "/" + local_bam_file,
    ), output_dir


@custom_task(cpu=24, memory=48, storage_gib=300)
def BamSort(
    bamfile: LatchFile, sort_by_name: bool, output_dir: LatchDir
) -> (LatchFile, LatchFile, LatchDir):
    # samtools sort -@ $nthread $inPath -o $outPath
    # samtools index -@ $nthread $outPath
    if sort_by_name == True:
        local_sorted_name = Path(bamfile).name.replace(".bam", ".sorted_read.bam")
        local_sorted_bamfile = Path(local_sorted_name).resolve()
        cmd_sort = [
            "samtools",
            "sort",
            "-@",
            "24",
            "-n",
            str(bamfile.local_path),
            "-o",
            str(local_sorted_bamfile),
        ]
        subprocess.run(" ".join(cmd_sort), shell=True, check=True)
        cmd_idx = ["touch", str(local_sorted_bamfile) + ".bai"]
        subprocess.run(" ".join(cmd_idx), shell=True, check=True)
    else:
        local_sorted_name = Path(bamfile).name.replace(".bam", ".sorted.bam")
        local_sorted_bamfile = Path(local_sorted_name).resolve()
        cmd_sort = [
            "samtools",
            "sort",
            "-@",
            "24",
            str(bamfile.local_path),
            "-o",
            str(local_sorted_bamfile),
        ]
        subprocess.run(" ".join(cmd_sort), shell=True, check=True)
        cmd_idx = ["samtools", "index", "-@", "24", str(local_sorted_bamfile)]
        subprocess.run(" ".join(cmd_idx), shell=True, check=True)
    print(" ".join(cmd_sort))

    return (
        LatchFile(
            local_sorted_bamfile, output_dir.remote_directory + "/" + local_sorted_name
        ),
        LatchFile(
            str(local_sorted_bamfile) + ".bai",
            output_dir.remote_directory + "/" + local_sorted_name + ".bai",
        ),
        output_dir,
    )


@custom_task(cpu=8, memory=16, storage_gib=100)
def Picard_RmDuplicates(
    input_bam: LatchFile, output_dir: LatchDir
) -> (LatchFile, LatchFile, LatchDir):
    bamfile_name = Path(input_bam.local_path).name.replace(".bam", "")
    out_bam = bamfile_name + ".nodups.bam"
    met_file = bamfile_name + ".nodups.metrics"

    cmd_picard = [
        "java",
        "-jar",
        "/root/picard.jar",
        "MarkDuplicates",
        "-I",
        input_bam.local_path,
        "-O",
        out_bam,
        "-M",
        met_file,
        "--REMOVE_DUPLICATES",
        "true",
    ]
    subprocess.run(" ".join(cmd_picard), shell=True, check=True)

    return (
        LatchFile(out_bam, output_dir.remote_directory + "/" + out_bam),
        LatchFile(met_file, output_dir.remote_directory + "/" + met_file),
        output_dir,
    )


@custom_task(cpu=4, memory=8, storage_gib=100)
def Shift_Tn5(
    bamfile: LatchFile, baifile: LatchFile, output_dir: LatchDir
) -> (LatchFile, LatchDir):
    copy(bamfile.local_path, str(Path().resolve()))
    copy(baifile.local_path, str(Path().resolve()))

    local_bamfile = Path(bamfile.local_path).name.replace(".bam", ".shifted.bam")
    cmd_shift = [
        "alignmentSieve",
        "-v",
        "-b",
        Path(bamfile.local_path).name,
        "-o",
        local_bamfile,
        "--ATACshift",
    ]
    subprocess.run(" ".join(cmd_shift), shell=True, check=True)
    return (
        LatchFile(local_bamfile, str(output_dir.remote_path) + "/" + local_bamfile),
        output_dir,
    )


@custom_task(cpu=4, memory=200, storage_gib=100)
def Run_Rscript_ATACSeqQC(
    bamfile: LatchFile,
    baifile: LatchFile,
    shifted_bamfile: LatchFile,
    shifted_baifile: LatchFile,
    saturation_bamfile: LatchFile,
    saturation_baifile: LatchFile,
    output_dir: LatchDir,
) -> LatchDir:
    local_dir = "BamQC/"
    copy(bamfile.local_path, str(Path().resolve()))
    copy(baifile.local_path, str(Path().resolve()))
    copy(shifted_bamfile.local_path, str(Path().resolve()))
    copy(shifted_baifile.local_path, str(Path().resolve()))
    copy(saturation_bamfile.local_path, str(Path().resolve()))
    copy(saturation_baifile.local_path, str(Path().resolve()))
    if not isdir(local_dir):
        os.mkdir(local_dir)
    bamfile.download()
    shifted_bamfile.download()
    saturation_bamfile.download()

    cmd_RunRscript = [
        "mamba",
        "run",
        "-n",
        "atacseqqc",
        "Rscript",
        "/root/wf/ATACSeqQC_Plots.r",
        Path(bamfile.local_path).name,
        Path(shifted_bamfile.local_path).name,
        Path(saturation_bamfile.local_path).name,
        local_dir,
    ]
    subprocess.run(" ".join(cmd_RunRscript), shell=True, check=True)
    return LatchDir(local_dir, output_dir.remote_directory + "/BamQC/")


@custom_task(cpu=4, memory=8, storage_gib=100)
def Call_Peaks(
    bamfile: LatchFile,
    output_dir: LatchDir,
    extSize: int = 150,
    shift: int = -75,
    slocal: int = 5000,
    llocal: int = 20000,
    p_val: float = 0.01,
) -> (LatchDir, LatchDir):
    base_name = Path(bamfile.local_path).name.replace(".sorted.bam", "")
    local_directory = base_name + "/MACS2/"
    cmd_macs2 = [
        "mamba",
        "run",
        "-n",
        "MACS2",
        "macs2",
        "callpeak",
        "-t",
        bamfile.local_path,
        "-g",
        "hs",
        "--nomodel",
        "--extsize",
        str(extSize),
        "--shift",
        str(shift),
        "--slocal",
        str(slocal),
        "--llocal",
        str(llocal),
        "-B",
        "--SPMR",
        "--keep-dup",
        "all",
        "-p",
        str(p_val),
        "-f",
        "BAMPE",
        "--outdir",
        local_directory,
        "--name",
        base_name,
    ]
    subprocess.run(" ".join(cmd_macs2), shell=True, check=True)
    return (
        LatchDir(local_directory, output_dir.remote_directory + "/MACS2/"),
        output_dir,
    )


@custom_task(cpu=8, memory=8, storage_gib=100)
def Bam2BigWig(
    bamfile: LatchFile, baifile: LatchFile, output_dir: LatchDir
) -> (LatchFile, LatchDir):
    copy(bamfile.local_path, str(Path().resolve()))
    copy(baifile.local_path, str(Path().resolve()))

    local_bamfile = Path(bamfile).name
    base_name = Path(bamfile).name.replace(".bam", "")
    local_file = base_name + ".bw"
    cmd_bigwig = ["bamCoverage", "-b", local_bamfile, "-o", local_file]
    subprocess.run(" ".join(cmd_bigwig), shell=True, check=True)
    return (
        LatchFile(local_file, output_dir.remote_directory + "/" + local_file),
        output_dir,
    )


@custom_task(cpu=8, memory=8, storage_gib=100)
def SamtoolsMerge(
    bamfile1: LatchFile, bamfile2: LatchFile, output_dir: LatchDir, sample: str
) -> (LatchFile, LatchFile, LatchDir):
    local_file = sample + ".merged.bam"
    cmd_merge = [
        "samtools",
        "merge",
        "-@",
        "8",
        local_file,
        bamfile1.local_path,
        bamfile2.local_path,
    ]
    subprocess.run(" ".join(cmd_merge), shell=True, check=True)

    cmd_idx = ["samtools", "index", "-@", "8", str(local_file)]
    print(" ".join(cmd_idx))
    subprocess.run(" ".join(cmd_idx), shell=True, check=True)
    print(listdir(Path()))
    return (
        LatchFile(local_file, output_dir.remote_directory + "/" + local_file),
        LatchFile(
            local_file + ".bai", output_dir.remote_directory + "/" + local_file + ".bai"
        ),
        output_dir,
    )


@custom_task(cpu=4, memory=4, storage_gib=100)
def Create_Pseudo_Replicates(
    output_dir: LatchDir,
    sample_id: str,
    fastq_fow_rep_1: LatchFile,
    fastq_rev_rep_1: LatchFile,
    fastq_fow_rep_2: Optional[LatchFile] = None,
    fastq_rev_rep_2: Optional[LatchFile] = None,
) -> (LatchFile, LatchFile, LatchFile, LatchFile):
    if (fastq_fow_rep_2 is None) or (fastq_rev_rep_2 is None):
        localdir = str(Path().resolve())
        cmd_split = [
            "seqkit",
            "split2",
            "-1",
            fastq_fow_rep_1.local_path,
            "-2",
            fastq_rev_rep_1.local_path,
            "-p",
            "2",
            "-O",
            localdir,
        ]
        subprocess.run(" ".join(cmd_split), shell=True, check=True)
        move(
            Path(fastq_fow_rep_1.local_path).name.replace(".fastq.gz", "")
            + ".part_001.fastq.gz",
            sample_id + ".pseudo_replicate_1.R1.fastq.gz",
        )
        move(
            Path(fastq_fow_rep_1.local_path).name.replace(".fastq.gz", "")
            + ".part_002.fastq.gz",
            sample_id + ".pseudo_replicate_2.R1.fastq.gz",
        )
        move(
            Path(fastq_rev_rep_1.local_path).name.replace(".fastq.gz", "")
            + ".part_001.fastq.gz",
            sample_id + ".pseudo_replicate_1.R2.fastq.gz",
        )
        move(
            Path(fastq_rev_rep_1.local_path).name.replace(".fastq.gz", "")
            + ".part_002.fastq.gz",
            sample_id + ".pseudo_replicate_2.R2.fastq.gz",
        )
        return (
            LatchFile(
                sample_id + ".pseudo_replicate_1.R1.fastq.gz",
                output_dir.remote_directory
                + "/"
                + sample_id
                + ".pseudo_replicate_1.R1.fastq.gz",
            ),
            LatchFile(
                sample_id + ".pseudo_replicate_2.R1.fastq.gz",
                output_dir.remote_directory
                + "/"
                + sample_id
                + ".pseudo_replicate_2.R1.fastq.gz",
            ),
            LatchFile(
                sample_id + ".pseudo_replicate_1.R2.fastq.gz",
                output_dir.remote_directory
                + "/"
                + sample_id
                + ".pseudo_replicate_1.R2.fastq.gz",
            ),
            LatchFile(
                sample_id + ".pseudo_replicate_2.R2.fastq.gz",
                output_dir.remote_directory
                + "/"
                + sample_id
                + ".pseudo_replicate_2.R2.fastq.gz",
            ),
        )
    else:
        return (fastq_fow_rep_1, fastq_fow_rep_2, fastq_rev_rep_1, fastq_rev_rep_2)


@custom_task(cpu=4, memory=8, storage_gib=100)
def Sort_Narrow_Peaks(MACS2_dir: LatchDir, tag: str) -> LatchFile:
    files = MACS2_dir.iterdir()
    for f in files:
        if Path(f).name.endswith("." + tag):
            narrow_peakfile = f
            break
    local_filename = Path(narrow_peakfile).name.replace("." + tag, ".sorted." + tag)
    cmd_sort_bed = ["sort-bed", narrow_peakfile.local_path, ">", local_filename]
    subprocess.run(" ".join(cmd_sort_bed), shell=True, check=True)
    return LatchFile(local_filename, MACS2_dir.remote_directory + "/" + local_filename)


@custom_task(cpu=4, memory=8, storage_gib=100)
def Run_IDR(
    sorted_peaks_rep1: LatchFile,
    sorted_peaks_rep2: LatchFile,
    oracle_peaks: LatchFile,
    sample: str,
    output_dir: LatchDir,
    idr_thresh: float = 0.01,
) -> (LatchFile, LatchFile, LatchFile):
    output_file = sample + ".idrValues"
    error_file = sample + "idrError.txt"
    cmd_idr = [
        "mamba",
        "run",
        "-n",
        "idr",
        "idr",
        "--verbose",
        "--samples",
        sorted_peaks_rep1.local_path,
        sorted_peaks_rep2.local_path,
        "--peak-list",
        oracle_peaks.local_path,
        "--rank",
        "p.value",
        "--plot",
        "--soft-idr-threshold",
        str(idr_thresh),
        "-o",
        output_file,
        "-l",
        error_file,
    ]
    subprocess.run(" ".join(cmd_idr), shell=True, check=True)
    cmd_sort_bed = ["sort-bed", output_file, ">", sample + ".sorted.idrValues"]
    subprocess.run(" ".join(cmd_sort_bed), shell=True, check=True)
    return (
        LatchFile(output_file, output_dir.remote_directory + "/" + output_file),
        LatchFile(error_file, output_dir.remote_directory + "/" + error_file),
        LatchFile(
            sample + ".sorted.idrValues",
            output_dir.remote_directory + "/" + sample + ".sorted.idrValues",
        ),
    )


@custom_task(cpu=4, memory=48, storage_gib=50)
def Compress_Coverages_Sample(
    bigWig_file: LatchFile, sample: str, outPath: LatchDir
) -> LatchDir:
    local_dir = "cov_parquet"
    if os.path.isdir(local_dir):
        os.mkdir(local_dir)

    print(bigWig_file, outPath, sample)
    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)

    wf.Compress_Coverages.Compress_Coverage(bigWig_file.local_path, local_dir, sample)

    return LatchDir(local_dir, os.path.join(outPath.remote_directory, local_dir))


### Sample: sample_id
### BamFile: merged_bam
### IDR_Peaks_Merged: idr_sorted+


@small_task
def Update_Registry_Tables(
    sample: str,
    bamfile: LatchFile,
    replicate_1_bamfile: LatchFile,
    replicate_2_bamfile: LatchFile,
    IDR_Peaks_Merged: LatchFile,
    run_name: str,
) -> str:
    target_project_name = "ATAC-Seq-Workflows"
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
    columns = [
        "sample",
        "bamfile",
        "replicate_1_bamfile",
        "replicate_2_bamfile",
        "IDR_Peaks_Merged",
    ]
    col_types = [str, LatchFile, LatchFile, LatchFile, LatchFile]
    values = [
        sample,
        bamfile,
        replicate_1_bamfile,
        replicate_2_bamfile,
        IDR_Peaks_Merged,
    ]

    d = dict(zip(columns, values))

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

    with target_table.update() as updater:
        ctr = len(target_table.get_dataframe())
        updater.upsert_record(name=str(ctr), **d)

    return target_table.id


@small_task
def Update_Registry_Tables_Plots(
    sample: str, output_dir: LatchDir, cov_path: LatchDir, run_name: str
) -> str:
    target_project_name = "ATAC-Seq-Workflows"
    target_table_name = f"{run_name}_Plots"
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
    columns = [
        "sample",
        "Fragment_Distribution",
        "Feature_Alignment_Coverage",
        "Saturation_Curves",
        "Peak_Coverage",
    ]
    col_types = [str, LatchFile, LatchFile, LatchFile, LatchDir]
    values = [
        sample,
        LatchFile(f"{output_dir.remote_directory}/Frag_Sizes.txt"),
        LatchFile(f"{output_dir.remote_directory}/featurealignment_coverage.txt"),
        LatchFile(f"{output_dir.remote_directory}/Saturation_Plots.txt"),
        cov_path,
    ]
    d = dict(zip(columns, values))

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
    with target_table.update() as updater:
        ctr = len(target_table.get_dataframe())
        updater.upsert_record(name=str(ctr), **d)

    return target_table.id

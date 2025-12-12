"""Module to retrieve read count statistics from log files."""

import os
from snakemake.script import snakemake
import pandas as pd


def parse_cutadapt_stats(cutadapt_logfile, is_paired):
    """Parse cutadapt log file to extract read statistics.

    Parameters
    ----------
    cutadapt_logfile : str
        Path to log file from ``cutadapt v5.1``
    is_paired : {True, False}
        ``True`` if reads are paired, ``False`` for single-end reads.

    Raises
    ------
    Exception
        If `cutadapt_logfile` is not properly formatted.

    Returns
    -------
    stats_dict : dict
        Keys:
            Total_raw_reads
            R1_reads_with_adapter
            R2_reads_with_adapter
            Total_trimmed_reads
    """
    with open(cutadapt_logfile, "r") as file:
        lines = file.readlines()

    # Check that log is properly formatted
    if "This is cutadapt 5.1" not in lines[0]:
        raise Exception(f"Error.. {cutadapt_logfile} is not properly formatted.")

    stats_dict = {
        "Total_raw_reads": None,
        "R1_reads_with_adapter": None,
        "R2_reads_with_adapter": None,
        "Total_trimmed_reads": None,
    }

    for line in lines:
        # Total reads
        if is_paired and "Total read pairs processed:" in line:
            stats_dict["Total_raw_reads"] = int(
                line.split(":")[1].strip().replace(",", "")
            )
        elif not is_paired and "Total reads processed:" in line:
            stats_dict["Total_raw_reads"] = int(
                line.split(":")[1].strip().replace(",", "")
            )
        # Reads with adapter
        if is_paired and "Read 1 with adapter:" in line:
            stats_dict["R1_reads_with_adapter"] = int(
                line.split(":")[1].split("(")[0].strip().replace(",", "")
            )
        elif not is_paired and "Reads with adapters:" in line:
            stats_dict["R1_reads_with_adapter"] = int(
                line.split(":")[1].split("(")[0].strip().replace(",", "")
            )
        if "Read 2 with adapter:" in line:
            stats_dict["R2_reads_with_adapter"] = int(
                line.split(":")[1].split("(")[0].strip().replace(",", "")
            )
        # Trimmed reads
        if is_paired and "Pairs written (passing filters):" in line:
            stats_dict["Total_trimmed_reads"] = int(
                line.split(":")[1].split("(")[0].strip().replace(",", "")
            )
        elif not is_paired and "Reads written (passing filters):" in line:
            stats_dict["Total_trimmed_reads"] = int(
                line.split(":")[1].split("(")[0].strip().replace(",", "")
            )

    if not is_paired:
        # No R2 reads
        stats_dict["R2_reads_with_adapter"] = 0

    # Check that all values were succesfully retrieved from log
    missing = [k for k, v in stats_dict.items() if v is None]
    if missing:
        raise Exception(f"Error.. Could not find {missing} in {cutadapt_logfile}.")

    return stats_dict


def parse_bbmerge_stats(bbmerge_logfile, trimmed_reads):
    """Parse bbmerge log file to extract read statistics.

    Parameters
    ----------
    bbmerge_logfile : str
        Path to log file from ``bbmerge-bbtools v39.52``
    trimmed_reads : int
        Number of trimmed reads

    Raises
    ------
    Exception
        If `bbmerge_logfile` is not properly formatted.
    ValueError
        If the number of reads in trimming output (`trimmed_reads`)
        does not correspond to the number of reads in merging input.

    Returns
    -------
    stats_dict : dict
        Keys:
            Total_merged_reads
            Nb_no_merging_solution
    """
    with open(bbmerge_logfile, "r") as file:
        lines = file.readlines()

    # Check that log is properly formatted
    if "Version 39.52" not in lines[-23]:
        raise Exception(f"Error.. {bbmerge_logfile} is not properly formatted.")

    # Check that number of processed reads corresponds to number of trimmed reads
    if "Pairs:" in lines[-16]:
        nb_reads = int(lines[-16].split()[1])
        if nb_reads != trimmed_reads:
            raise ValueError(
                "Error.. Number of written reads in trim output does not correspond to number of processed reads in merge input"
            )

    stats_dict = {"Total_merged_reads": None, "Nb_no_merging_solution": None}

    if "Joined:" in lines[-15]:
        stats_dict["Total_merged_reads"] = int(lines[-15].split()[1])

    if "No Solution:" in lines[-13]:
        stats_dict["Nb_no_merging_solution"] = int(
            lines[-13].split("No Solution")[1].split()[1]
        )

    # Check that all values were succesfully retrieved from log
    missing = [k for k, v in stats_dict.items() if v is None]
    if missing:
        raise Exception(f"Error.. Could not find {missing} in {bbmerge_logfile}.")

    return stats_dict


def parse_vsearch_stats(vsearch_logfile, merged_reads):
    """Parse vsearch log file to extract read statistics.

    Parameters
    ----------
    vsearch_logfile : str
        Path to log file from ``vsearch v2.29.3``
    merged_reads : int
        Number of merged reads

    Raises
    ------
    Exception
        If `vsearch_logfile` is not properly formatted.
    ValueError
        If the number of reads in merging output (`merged_reads`)
        does not correspond to the number of reads in aggregating input.

    Returns
    -------
    stats_dict : dict
        Keys:
            Total_aggregated_reads
            Nb_singletons
    """
    with open(vsearch_logfile, "r") as file:
        lines = file.readlines()

    # Check that log is properly formatted
    if "vsearch v2.29.3" not in lines[0]:
        raise Exception(f"Error.. {vsearch_logfile} is not properly formatted.")

    # Check that number of processed reads corresponds to number of merged reads
    if " seqs, min " in lines[-5]:
        nb_reads = int(lines[-5].split(" seqs, min ")[0].split(" nt in ")[1])
        if nb_reads != merged_reads:
            raise ValueError(
                "Error.. Number of written reads in merge output does not correspond to number of processed reads in aggregate input"
            )

    stats_dict = {}

    # Add singletons total
    if "clusters discarded" in lines[-1]:
        singletons = lines[-1].split(" clusters discarded")[0].split(",")[-1].strip()
    else:
        singletons = 0
    stats_dict["Nb_singletons"] = int(singletons)

    # Calculate number of reads corresponding to aggregated sequences
    stats_dict["Total_aggregated_reads"] = merged_reads - stats_dict["Nb_singletons"]

    return stats_dict


def generate_read_stats(
    cutadapt_logfile,
    bbmerge_logfile,
    vsearch_logfile,
    Nout_logfile,
    outpath,
    sample_name,
    is_paired,
):
    r"""Parse different types of log files and extract read count statistics.

    Parameters
    ----------
    cutadapt_logfile : str
        Path to log file from ``cutadapt v5.1``
    bbmerge_logfile : str
        Path to log file from ``bbmerge-bbtools v39.52``
    vsearch_logfile : str
        Path to log file from ``vsearch v2.29.3``
    Nout_logfile : str
        Path to text file containing the number of DNA sequencing reads
        discarded because they contained N(s).
    outpath : str
        Path to save output dataframe of read count statistics.
    sample_name : str
        Sample identifier
    is_paired : {True, False}
        ``True`` if reads are paired, ``False`` for single-end reads.

    Notes
    -----
    Statistics retrieved include:

    * the total number of reads
    * the number of trimmed reads
    * the number of merged reads
    * the number of aggregated reads
    """
    # Step 1 - Parse cutadapt log
    cutadapt_stats = parse_cutadapt_stats(cutadapt_logfile, is_paired)
    fullstats = pd.DataFrame({sample_name: cutadapt_stats}).T

    # Step 2 - Parse bbmerge log if not empty
    if os.path.getsize(bbmerge_logfile) > 0:
        bbmerge_stats = parse_bbmerge_stats(
            bbmerge_logfile, trimmed_reads=cutadapt_stats["Total_trimmed_reads"]
        )
        for k, v in bbmerge_stats.items():
            fullstats[k] = v
    else:
        # case when single-end
        fullstats["Total_merged_reads"] = cutadapt_stats["Total_trimmed_reads"]

    # Step 3 - Parse vsearch log
    vsearch_stats = parse_vsearch_stats(
        vsearch_logfile,
        merged_reads=fullstats.loc[sample_name, "Total_merged_reads"],
    )
    for k, v in vsearch_stats.items():
        fullstats[k] = v

    # Step 4 - Add number of discarded sequences (with Ns) at FASTA parsing step
    Nout_df = pd.read_csv(Nout_logfile, header=None)
    fullstats["Contain_Ns"] = Nout_df.iat[0, 1]

    ### Derive metrics

    # Calculate percentage of trimmed reads
    fullstats["Trimmed_%"] = (
        fullstats["Total_trimmed_reads"] / fullstats["Total_raw_reads"]
    )

    # Calculate percentage of properly merged reads relative to number of trimmed reads
    fullstats["Merged_%"] = (
        fullstats["Total_merged_reads"] / fullstats["Total_trimmed_reads"]
    )

    # Calculate percentage of reads corresponding to non-singletons (relative to number of properly merged reads)
    fullstats["Aggregated_%"] = (
        fullstats["Total_aggregated_reads"] / fullstats["Total_merged_reads"]
    )

    ### Output to csv file
    fullstats[
        [
            "Total_raw_reads",
            "R1_reads_with_adapter",
            "R2_reads_with_adapter",
            "Total_trimmed_reads",
            "Trimmed_%",
            "Nb_no_merging_solution",
            "Total_merged_reads",
            "Merged_%",
            "Nb_singletons",
            "Total_aggregated_reads",
            "Aggregated_%",
            "Contain_Ns",
        ]
    ].reset_index(names="Sample_name").to_csv(outpath, index=False)

    return


generate_read_stats(
    snakemake.input.cutadapt_log,
    snakemake.input.bbmerge_log,
    snakemake.input.vsearch_log,
    snakemake.input.N_discarded_log,
    snakemake.output[0],
    snakemake.wildcards.sample,
    snakemake.params.is_paired,
)

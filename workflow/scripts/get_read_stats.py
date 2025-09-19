"""Module to retrieve read count statistics from log files."""

import os
from snakemake.script import snakemake
import pandas as pd


def generate_read_stats(
    cutadapt_logfile,
    pandaseq_logfile,
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
    pandaseq_logfile : str
        Path to log file from ``pandaseq v2.11``
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
    
    Raises
    ------
    Exception
        If any log file is not properly formatted.
    
    Notes
    -----
    Statistics retrieved include:

    * the total number of reads
    * the number of trimmed reads
    * the number of merged reads
    * the number of aggregated reads
    * the number of reads lost at each step
    """
    stats_dict = {}

    # Step 1 - Process cutadapt log

    # Parse cutadapt stats
    with open(cutadapt_logfile, "r") as file:
        lines = file.readlines()

    # Check that log is properly formatted
    if "This is cutadapt 5.1" not in lines[0]:
        raise Exception(
            f"Error.. {cutadapt_logfile} is not properly formatted. Double-check cutadapt's wrapper version."
        )

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

    # Create dataframe to store all stats
    fullstats = pd.DataFrame(
        {sample_name: stats_dict},
        index=[
            "Total_raw_reads",
            "R1_reads_with_adapter",
            "R2_reads_with_adapter",
            "Total_trimmed_reads",
        ],
    ).T  # transpose

    # Calculate percentage of trimmed reads
    fullstats["Trimmed_%"] = (
        fullstats["Total_trimmed_reads"] / fullstats["Total_raw_reads"]
    )

    # Calculate number of reads discarded at the trimming stage
    fullstats["Trimming"] = (
        fullstats["Total_raw_reads"] - fullstats["Total_trimmed_reads"]
    )

    # Step 2 - Process pandaseq log

    # Check that log is not empty
    if os.path.getsize(pandaseq_logfile) > 0:
        with open(pandaseq_logfile, "r") as file:
            first_line = file.readline()

        # Check log header
        if "INFO	VER	pandaseq 2.11" not in first_line:
            raise Exception(
                f"Error.. {pandaseq_logfile} is not properly formatted. Make sure you've added the --use-conda flag in the snakemake command line, which specifies the correct package versions to be used"
            )

        # Parse pandaseq stats
        logfile = pd.read_csv(
            pandaseq_logfile,
            sep="\t",
            skiprows=21,
            skipfooter=1,
            engine="python",
            names=["id", "err_stat", "field", "value", "details"],
        )
        stats = logfile[
            logfile.field.isin(["LOWQ", "NOALGN", "OK", "READS", "SLOW"])
        ].iloc[-5:, :][["field", "value"]]
        stats["value"] = stats.value.astype(int)

        # Validation step - check that the number of processed reads corresponds from trim output to merge input
        if (
            fullstats.loc[fullstats.index == sample_name, "Total_trimmed_reads"].item()
            != stats.loc[stats.field == "READS", "value"].item()
        ):
            print(
                "---ERROR---\nNumber of written reads in trim output does not correspond to number of processed reads in merge input"
            )

        # Add stats to dataframe
        fullstats.loc[
            fullstats.index == sample_name,
            ["Not_merged_LOWQ", "Not_merged_NOALGN", "Total_merged_reads"],
        ] = (
            stats.set_index("field").T[["LOWQ", "NOALGN", "OK"]].values
        )

        # Cast type to remove ','
        fullstats = fullstats.astype(
            {
                "Not_merged_LOWQ": "int",
                "Not_merged_NOALGN": "int",
                "Total_merged_reads": "int",
            }
        )

    else:  # Empty log likely touched when trimming single-end reads = no need for merging with pandaseq
        fullstats["Total_merged_reads"] = fullstats["Total_trimmed_reads"]

    # Calculate percentage of properly merged reads relative to number of trimmed reads
    fullstats["Merged_%"] = (
        fullstats["Total_merged_reads"] / fullstats["Total_trimmed_reads"]
    )

    # Calculate number of reads discarded at the merging step
    fullstats["Merging"] = (
        fullstats["Total_trimmed_reads"] - fullstats["Total_merged_reads"]
    )

    # Step 3 - Process vsearch log

    # Parse vsearch stats
    with open(vsearch_logfile, "r") as file:
        lines = file.readlines()

    # Check that log is properly formatted
    if "vsearch v2.29.3" not in lines[0]:
        raise Exception(
            f"Error.. {vsearch_logfile} is not properly formatted. Make sure you've added the --use-conda flag in the snakemake command line, which specifies the correct package versions to be used"
        )

    # Add singletons total
    if "clusters discarded" in lines[-1]:
        singletons = lines[-1].split(" clusters discarded")[0].split(",")[-1].strip()
    else:
        singletons = 0
    fullstats.loc[fullstats.index == sample_name, "Nb_singletons"] = singletons

    # Convert column type to int
    fullstats["Nb_singletons"] = fullstats["Nb_singletons"].astype(int)

    # Calculate number of reads corresponding to aggregated sequences
    fullstats["Nb_non-singletons"] = (
        fullstats["Total_merged_reads"] - fullstats["Nb_singletons"]
    )

    # Calculate percentage of reads corresponding to non-singletons (relative to number of properly merged reads)
    fullstats["Aggregated_%"] = (
        fullstats["Nb_non-singletons"] / fullstats["Total_merged_reads"]
    )

    # Calculate number of reads discarded at the aggregating step (= number of singletons)
    fullstats["Aggregating"] = fullstats["Nb_singletons"]

    # Step 4 - Add number of discarded sequences (with Ns) at parsing step

    Nout_df = pd.read_csv(Nout_logfile, header=None)
    fullstats["Contain_Ns"] = Nout_df.iat[0, 1]

    # Output to csv file
    fullstats.reset_index(names="Sample_name").to_csv(outpath, index=False)

    return


generate_read_stats(
    snakemake.input.cutadapt_log,
    snakemake.input.pandaseq_log,
    snakemake.input.vsearch_log,
    snakemake.input.N_discarded_log,
    snakemake.output[0],
    snakemake.wildcards.sample,
    snakemake.params.is_paired,
)

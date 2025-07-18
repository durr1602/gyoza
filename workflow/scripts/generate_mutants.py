from snakemake.script import snakemake

# from scripts.my_functions import load_codon_dic, get_single_double, get_nt_seq
import pandas as pd
import numpy as np
import itertools


def load_codon_dic(table):
    """
    Small function to load the codon table and return a dictionary
    """
    codon_table = pd.read_csv(table, header=0)
    codon_table["codon"] = codon_table["codon"].str.upper()
    codon_dic = dict(zip(codon_table["codon"], codon_table["aminoacid"]))
    return codon_dic


def get_alt_codons(seq, codon_dic, mode="NNN"):
    """
    Based on a DNA sequence, the function returns two lists:
    1) A list containing all 0-based amino acid positions for the sequence
    2) A list containing all possible alternative codons (other than WT codon) at the matching positions
    For list 2, the mode defines which codons are acceptable: NNN by default, or NNK
    Codons are fetched in the provided codon table (dictionary)
    """

    if mode == "NNN":
        alt = [x for x in codon_dic.keys()]
    elif mode == "NNK":
        alt = [x for x in codon_dic.keys() if x[2] in ["G", "T"]]
    else:
        print("Pleae specify a correct mode: either NNN or NNK")

    pos_l = []
    var_l = []

    for i in range(0, len(seq), 3):
        list_var = [x for x in alt if x != seq[i : i + 3]]
        pos_l.append(i // 3)  # 0-based position (aa)
        var_l.append(list_var)  # list of possible codons other than WT

    return pos_l, var_l


def get_single_double(df, codon_dic, codon_mode):
    """
    Processes a minimal dataframe containing the wild-type nucleotide sequence to mutate.
    From the codon_mode, generates all single mutants and optionally double mutants.
    Returns dataframe in long format (one row per mutated codon).
    """
    if codon_mode not in ["NNN", "NNK", "NNN x NNN", "NNK x NNK"]:
        raise Exception(f"Error.. check the codon mode specified in the config.")
    elif codon_mode in ["NNN", "NNK"]:
        single_codon_mode = codon_mode
    else:
        single_codon_mode = codon_mode.split(" x ")[0]

    singles_compact = df.copy()
    singles_compact["pos"], singles_compact["alt_codons"] = zip(
        *singles_compact.WT_seq.apply(
            lambda x: get_alt_codons(x, codon_dic, single_codon_mode)
        )
    )

    singles_compact = singles_compact.explode(["pos", "alt_codons"])
    singles_df = singles_compact.explode("alt_codons")
    singles_df["mutations"] = singles_df.apply(
        lambda row: {row[f"pos"]: row[f"alt_codons"]}, axis=1
    )

    if codon_mode in ["NNN x NNN", "NNK x NNK"]:
        pairwise_df = df.copy()

        # Get pairwise combinations of codons that will be mutated
        pairwise_df["combination"] = pairwise_df.WT_seq.apply(
            lambda x: [
                x
                for x in list(itertools.combinations(range(0, len(x) // 3), 2))
                for _ in range(2)
            ]
        )
        doubles_compact_df = pairwise_df.explode("combination")
        doubles_compact_df["mutated_codon"] = np.tile(
            [1, 2], len(doubles_compact_df) // 2 + 1
        )[: len(doubles_compact_df)]
        doubles_compact_df["pos"] = doubles_compact_df.apply(
            lambda row: row.combination[row.mutated_codon - 1], axis=1
        )

        # Build a lookup dictionary for alternative codons at each position
        # We use the previously built dataframe
        alt_lookup = singles_df.groupby("pos")["alt_codons"].unique().to_dict()

        # Leverage dictionary to map alternative codons at each position
        doubles_compact_df["alt_codons"] = doubles_compact_df["pos"].map(alt_lookup)

        # Dataframe is pivoted only to be able to use pd.explode(), then later on melted to go back to long format
        doubles_piv = doubles_compact_df.pivot_table(
            index=["Mutated_seq", "WT_seq", "combination"],
            columns="mutated_codon",
            values=["pos", "alt_codons"],
            aggfunc="first",
        ).reset_index()
        doubles_piv.columns = [x[0] for x in doubles_piv.columns[:-4]] + [
            f"{x[0]}{x[1]}" for x in doubles_piv.columns[-4:]
        ]

        # Reshape
        doubles_exp1 = doubles_piv.explode("alt_codons1")
        doubles_exp2 = doubles_exp1.explode("alt_codons2")
        doubles_df = doubles_exp2.reset_index(drop=True)

        # Create dictionary of mutations, here we go with numpy to be as efficient as possible, even with very large datasets
        keys = np.stack([doubles_df["pos1"].values, doubles_df["pos2"].values], axis=1)
        vals = np.stack(
            [doubles_df["alt_codons1"].values, doubles_df["alt_codons2"].values], axis=1
        )
        doubles_df["mutations"] = [dict(zip(k, v)) for k, v in zip(keys, vals)]
        doubles_df.drop(
            ["pos1", "pos2", "alt_codons1", "alt_codons2", "combination"],
            axis=1,
            inplace=True,
        )

    else:
        doubles_df = pd.DataFrame()

    # Concatenate dataframes for single mutants (and double mutants if there are any)
    mutants_df = pd.concat(
        [singles_df.drop(["pos", "alt_codons"], axis=1), doubles_df], ignore_index=True
    )

    return mutants_df


def get_nt_seq(seq, mut_dic):
    """
    Reconstitutes a nucleotide sequence based on a dictionary of mutations,
    containing positions and alternative codons.
    """
    list_codons = [
        seq[i : i + 3]
        for i in range(0, len(seq), 3)  # Convert nucleotide sequence to list of codons
    ]

    for pos, mut_codon in mut_dic.items():
        if 0 <= pos < len(list_codons):
            list_codons[pos] = mut_codon

    return "".join(list_codons)


def generate_mutants(wtseq_path, outpath, mutated_seq, codon_table, codon_mode):
    """
    Generates all single mutants and optionally all double mutants,
    based on the provided wild-type sequence
    """
    # Create directory if does not already exist
    # outdir.mkdir(parents=True, exist_ok=True)

    # Load codon dictionary
    codon_dic = load_codon_dic(codon_table)

    # Load dataframe with all wild-type sequences
    all_wt = pd.read_csv(wtseq_path)
    wt_df = all_wt[all_wt["Mutated_seq"] == mutated_seq].copy()
    wt_df["WT_seq"] = wt_df["WT_seq"].str.upper()
    wt_df["nt_seq"] = wt_df["WT_seq"]

    # Get all single mutants (and double mutants if necessary)
    mutants_df = get_single_double(wt_df, codon_dic, codon_mode)

    # Reconstruct nucleotide sequences from mutations
    mutants_df["nt_seq"] = list(
        itertools.starmap(
            get_nt_seq, zip(mutants_df["WT_seq"], mutants_df["mutations"])
        )
    )
    # Concatenate with WT
    # The mutations are dropped so we can annotate anything in downstream steps
    expmut_df = pd.concat(
        [wt_df, mutants_df.drop(["mutations"], axis=1)], ignore_index=True
    )

    # Convert to slightly lighter data types
    for x in expmut_df.columns:
        expmut_df[x] = expmut_df[x].astype("string")

    # Export dataframe
    expmut_df.to_csv(outpath, index=False, compression="gzip", chunksize=100_000)

    return


generate_mutants(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.wildcards.mutseq,
    snakemake.config["codon"]["table"],
    snakemake.config["codon"]["mode"],
)

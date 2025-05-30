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
    *singles_compact.WT_seq.apply(lambda x: get_alt_codons(x, codon_dic, single_codon_mode))
    )

    singles_compact = singles_compact.explode(["pos", "alt_codons"])
    singles_df = singles_compact.explode("alt_codons")
    singles_df["mutations"] = singles_df.apply(lambda row: {row[f"pos"]: row[f"alt_codons"]}, axis=1)

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
        alt_lookup = (
            singles_df
            .groupby("pos")["alt_codons"]
            .unique()
            .to_dict()
        )

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
        doubles_df.drop(["pos1", "pos2", "alt_codons1", "alt_codons2", "combination"], axis=1, inplace=True)
    
    else:
        doubles_df = pd.DataFrame()

    # Concatenate dataframes for single mutants (and double mutants if there are any)
    mutants_df = pd.concat([singles_df.drop(["pos", "alt_codons"], axis=1), doubles_df], ignore_index=True)

    return mutants_df

def get_mutations(seq, wt, codon_dic):
    """
    By comparing a mutated DNA sequence to the wild-type sequence,
    this function returns the mutations (if there are any).
    Mutations are formatted as # mutated codon / position / alternative codon / alternative amino acid
    in lists with matching indexes to be able to quickly convert to 1 row per mutation per mutated codon
    The alternative and corresponding wild-type codons are translated into their corresponding amino acid using the provided codon table dictionary
    From there, we also calculate the Hamming distances in codons, nucleotides and amino acids.
    """
    if len(seq) != len(wt):
        raise ValueError(
            f"Error.. Cannot annotate expected mutants because at least one sequence is of different length than wild-type."
        )

    if len(seq) % 3 != 0:
        raise ValueError(
            f"Error.. the length of the DNA sequence is not a multiple of 3."
        )
    
    VALID_BASES = {"A", "C", "G", "T"}

    if not set(seq).issubset(VALID_BASES):
        raise ValueError(
        f"Error.. one of the provided nucleotide sequences contains an unrecognized character."
        )

    is_wt = seq == wt

    list_alt_pos, list_alt_cod, list_alt_aa, list_aa = [], [], [], []
    Nham_nt = 0
    Nham_aa = 0

    wt_codons = [
        wt[i : i + 3]
        for i in range(
            0, len(wt), 3
        )  # Converting WT nucleotide sequence to list of codons
    ]
    seq_codons = [
        seq[i : i + 3]
        for i in range(
            0, len(seq), 3
        )  # Converting nucleotide sequence of variant to list of codons
    ]

    for i, (wtc, c) in enumerate(zip(wt_codons, seq_codons)):  # Loop through codons
        alt_aa = codon_dic.get(c)
        wt_aa = codon_dic.get(wtc)
        list_aa.append(alt_aa)
        if c != wtc:
            list_alt_pos.append(i)
            list_alt_cod.append(c)
            list_alt_aa.append(alt_aa)
            Nham_nt += (wtc != c) * sum(
                x != y for x, y in zip(wtc, c)  # Calls zip only when codons differ
            )
            if alt_aa != wt_aa:
                Nham_aa += 1

    Nham_codons = len(list_alt_pos)

    if Nham_codons > 0:
        mut_codons = list(range(1, Nham_codons + 1))
    else:
        mut_codons = [0]
        list_alt_pos = ["not-applicable"]
        list_alt_cod = ["not-applicable"]
        list_alt_aa = ["not-applicable"]
    
    aa_seq = "".join(list_aa)

    return (
        is_wt,
        aa_seq,
        Nham_codons,
        Nham_nt,
        Nham_aa,
        mut_codons,
        list_alt_pos,
        list_alt_cod,
        list_alt_aa,
    )

def annotate_mutants(df, codon_dic):
    """
    Processes subdataframes with at least two columns (nt_seq and WT_seq)
    Compares both and returns mutations with a custom function
    Outputs corresponding subdataframe with annotated mutations
    """
    per_seq_cols = ["WT", "aa_seq", "Nham_codons", "Nham_nt", "Nham_aa"]
    per_mut_cols = ["mutated_codon", "pos", "alt_codons", "alt_aa"]
    new_cols = per_seq_cols + per_mut_cols

    # Making sure sequences are capitalized
    df["nt_seq"] = df["nt_seq"].str.upper()
    df["WT_seq"] = df["WT_seq"].str.upper()

    collected_mutations = [
        get_mutations(seq, wt, codon_dic)
        for seq, wt in zip(df["nt_seq"], df["WT_seq"])
    ]

    mutations_dict = dict(zip(new_cols, zip(*collected_mutations)))
    df = df.assign(**mutations_dict)
    df = df.explode(per_mut_cols).reset_index(drop=True)

    return df


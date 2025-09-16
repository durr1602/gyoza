"""Module to annotate mutants.

Compares all mutated DNA sequences of a sample to the wild-type sequence.

Sequences must be the same length as the WT.
Mutations are collected by comparing codons at matching positions.

A codon table (genetic code) is interrogated to collect the corresponding
protein-level mutations.

"""

from snakemake.script import snakemake
import pandas as pd
import json


def load_codon_dic(table):
    r"""Convert the CSV-formatted codon table to a dict.
    
    Parameters
    ----------
    table : str
        Path to CSV-formatted codon table.
        Header must be on first line and include columns "codon" and "aminoacid"
    
    Returns
    -------
    dict
    """
    codon_table = pd.read_csv(table, header=0)
    codon_table["codon"] = codon_table["codon"].str.upper()
    codon_dic = dict(zip(codon_table["codon"], codon_table["aminoacid"]))
    return codon_dic


def get_mutations(seq, wt, codon_dic):
    r"""Collect differences between a mutated DNA sequence and the wild-type.
    
    Parameters
    ----------
    seq : str
        DNA sequence of the mutant (with bases either A, C, G or T).
        Length should be the same as `wt` and be a multiple of 3.
    wt : str
        Wild-type DNA sequence (with bases either A, C, G or T).
        Length should be the same as `seq` and be a multiple of 3.
    codon_dic : dict
        Codon table associating codons to amino acid residues.
    
    Returns
    -------
    is_wt : bool
        True if `seq` is `wt`
    aa_seq : str
        Protein sequence translated from `seq`
    Nham_codons : int
        Number of codon changes
    Nham_nt : int
        Number of nucleotide changes
    Nham_aa : int
        Number of amino acid changes
    protein_pos : list
        Positions in the protein sequence of each non-silent mutation
    protein_alt_aa : list
        Alternative residue in the protein sequence for each non-silent mutation
    mutated_codon : list
        1-based indexes for each mutated codon (nth mutated codon)
    mutation_pos : list
        Positions in `seq` of each mutation
    mutation_alt_codon : list
        Alternative codon in `seq` for each mutation
    mutation_alt_aa : list
        Alternative residue translated from `mutation_alt_codon`
    mutation_type : list
        Either "silent", "missense" or "nonsense" based on `mutation_alt_aa`
    
    Notes
    --------
    Mutations are formatted as # mutated codon / position / alternative codon /
    alternative amino acid, in lists with matching indexes to be able
    to quickly convert to 1 row per mutation per mutated codon.

    The alternative and corresponding wild-type codons are translated into
    their corresponding amino acid using the `codon_dic`.

    Sequence-level attributes include the Hamming distances (Nham),
    i.e. the number of codon, nucleotide and amino acid changes.
    """
    if len(seq) % 3 != 0:
        raise ValueError(
            f"Error.. the length of the DNA sequence is not a multiple of 3."
        )

    # Note: the two following checks are validated early on
    # but let's keep them just in case

    if len(seq) != len(wt):
        raise ValueError(
            f"Error.. Cannot annotate expected mutants because at least one sequence is of different length than wild-type."
        )

    if not set(seq).issubset({"A", "C", "G", "T"}):
        raise ValueError(
            f"Error.. one of the provided nucleotide sequences contains an unrecognized character."
        )

    is_wt = seq == wt

    mutation_pos, mutation_alt_codon, mutation_alt_aa, mutation_type = [], [], [], []
    full_aa_seq, full_wt_aa = [], []
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
        wt_aa = codon_dic.get(wtc)
        alt_aa = codon_dic.get(c)
        full_aa_seq.append(alt_aa)
        full_wt_aa.append(wt_aa)

        if c != wtc:
            mutation_pos.append(i)
            mutation_alt_codon.append(c)
            mutation_alt_aa.append(alt_aa)
            Nham_nt += sum(x != y for x, y in zip(wtc, c))
            if alt_aa != wt_aa:
                Nham_aa += 1
                mutation_type.append("nonsense" if alt_aa == "*" else "missense")
            else:
                mutation_type.append("silent")

    # Protein-level comparison
    aa_seq = "".join(full_aa_seq)
    wt_aa_seq = "".join(full_wt_aa)

    protein_pos, protein_alt_aa = [], []

    for i, (wt_aa, alt_aa) in enumerate(zip(wt_aa_seq, aa_seq)):
        if wt_aa != alt_aa:
            protein_pos.append(int(i)),
            protein_alt_aa.append(alt_aa)

    Nham_codons = len(mutation_pos)

    if Nham_codons > 0:
        mutated_codon = list(
            range(1, Nham_codons + 1)
        )  # 1-based index of mutated codons
        if Nham_aa == 0:
            protein_alt_aa = ["not-applicable"]
    else:  # Handle WT (no mutations)
        mutated_codon = [0]
        mutation_pos = ["not-applicable"]
        mutation_alt_codon = ["not-applicable"]
        mutation_alt_aa = ["not-applicable"]
        mutation_type = ["wt"]
        protein_alt_aa = ["not-applicable"]

    return (
        is_wt,
        aa_seq,
        Nham_codons,
        Nham_nt,
        Nham_aa,
        protein_pos,
        protein_alt_aa,
        mutated_codon,
        mutation_pos,
        mutation_alt_codon,
        mutation_alt_aa,
        mutation_type,
    )


def annotate_mutants(df, codon_dic):
    r"""Annotate a dataframe of mutated DNA sequences with mutations.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns: "nt_seq" and "WT_seq"
    codon_dic : dict
        Codon table associating codons to amino acid residues.
    
    Returns
    -------
    pandas.DataFrame
        `df` with additional columns to describe mutations.
    
    Notes
    ----------
    Uses custom function to collect mutations by comparing each sequence
    to the corresponding wild-type on the matching row.
    """
    per_seq_cols = [
        "WT",
        "aa_seq",
        "Nham_codons",
        "Nham_nt",
        "Nham_aa",
        "pos",
        "alt_aa",
    ]
    per_mut_cols = [
        "mutated_codon",
        "mutation_pos",
        "mutation_alt_codons",
        "mutation_alt_aa",
        "mutation_type",
    ]
    new_cols = per_seq_cols + per_mut_cols

    # Making sure sequences are capitalized
    df["nt_seq"] = df["nt_seq"].str.upper()
    df["WT_seq"] = df["WT_seq"].str.upper()

    collected_mutations = [
        get_mutations(seq, wt, codon_dic) for seq, wt in zip(df["nt_seq"], df["WT_seq"])
    ]

    mutations_dict = dict(zip(new_cols, zip(*collected_mutations)))
    df = df.assign(**mutations_dict)
    df = df.explode(per_mut_cols).reset_index(drop=True)

    return df


def get_annotated_mutants(mut_path, outpath, position_offset, codon_table):
    r"""Annotates non-empty dataframes of mutated DNA sequences (1 df per sample).

    Parameters
    ----------
    mut_path : str
        Path to dataframe containing mutated DNA sequences (and wild-type).
        A single wild-type sequence (repeated on every row) is expected.
    outpath : str
        Path to save output dataframe, either with additional columns or
        with column headers only (if input df is empty).
    position_offset : int
        Starting position in the full protein sequence.
    codon_table : str
        Path to CSV-formatted codon table.
        Header must be on first line and include columns "codon" and "aminoacid"
    
    Notes
    -------
    One column, "aa_pos" is obtained by adding the position offset to the
    positions of every non-silent mutation at the protein level.
    This column is left "unexploded", such that each row contains a list.
    Special care is taken to ensure that these lists can be parsed later on
    even if they are saved as strings in the CSV file at `outpath`.
    """
    # Load codon dictionary
    codon_dic = load_codon_dic(codon_table)

    # Load dataframe
    df = pd.read_csv(mut_path)

    # Annotate non-indels if there are any
    if not df.empty:
        annot_df = annotate_mutants(df, codon_dic)
        # Add position offset at mutation level
        annot_df["mutation_aa_pos"] = annot_df.mutation_pos.apply(
            lambda x: (
                int(x) + position_offset if x != "not-applicable" else "not-applicable"
            )
        )
        # Add position offset for amino acid changes at protein level
        # We also need to robustly convert the type and serialize
        # for proper parsing downstream
        annot_df["aa_pos"] = annot_df.pos.apply(
            lambda x: (
                json.dumps([int(y + position_offset) for y in x])
                if x != []
                else json.dumps(["not-applicable"])
            )
        )
        annot_df["alt_aa"] = annot_df.alt_aa.apply(json.dumps)
        # Remove obsolete columns
        annot_df.drop(["mutation_pos", "pos"], axis=1, inplace=True)
    else:
        # Rescue expected column headers
        annot_df = pd.DataFrame(
            columns=[
                "nt_seq",
                "WT_seq",
                "WT",
                "aa_seq",
                "Nham_codons",
                "Nham_nt",
                "Nham_aa",
                "aa_pos",
                "alt_aa",
                "mutated_codon",
                "mutation_aa_pos",
                "mutation_alt_codons",
                "mutation_alt_aa",
                "mutation_type",
            ]
        )

    # Save outputs
    annot_df.sort_values(by="WT", ascending=False).to_csv(outpath, index=False)

    return


get_annotated_mutants(
    snakemake.input[0],
    snakemake.output.annot_rc,
    snakemake.params.position_offset,
    snakemake.params.genetic_code,
)

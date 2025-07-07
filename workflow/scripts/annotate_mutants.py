from snakemake.script import snakemake

# from scripts.my_functions import load_codon_dic, annotate_mutants
import pandas as pd


def load_codon_dic(table):
    """
    Small function to load the codon table and return a dictionary
    """
    codon_table = pd.read_csv(table, header=0)
    codon_table["codon"] = codon_table["codon"].str.upper()
    codon_dic = dict(zip(codon_table["codon"], codon_table["aminoacid"]))
    return codon_dic


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
        get_mutations(seq, wt, codon_dic) for seq, wt in zip(df["nt_seq"], df["WT_seq"])
    ]

    mutations_dict = dict(zip(new_cols, zip(*collected_mutations)))
    df = df.assign(**mutations_dict)
    df = df.explode(per_mut_cols).reset_index(drop=True)

    return df


def get_annotated_mutants(mut_path, outpath, indel_outpath, codon_table):
    """
    Annotates mutants (input = 1 dataframe per sample).
    """

    # Load codon dictionary
    codon_dic = load_codon_dic(codon_table)

    # Load dataframe
    df = pd.read_csv(mut_path)

    # Identify indels
    is_indel = df["nt_seq"].str.len() != df["WT_seq"].str.len()
    df_indels = df[is_indel].copy()
    df_valid = df[~is_indel].copy()

    # Annotate non-indels if there are any
    if not df_valid.empty:
        annot_df = annotate_mutants(df_valid, codon_dic)
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
                "mutated_codon",
                "pos",
                "alt_codons",
                "alt_aa",
            ]
        )

    # Ensure dataframe with indels is created even if empty
    if df_indels.empty:
        df_indels = pd.DataFrame(columns=df.columns)

    # Save outputs
    annot_df.to_csv(outpath, index=False)
    df_indels.to_csv(indel_outpath, index=False)

    return


get_annotated_mutants(
    snakemake.input[0],
    snakemake.output.annot_rc,
    snakemake.output.indels,
    snakemake.config["codon"]["table"],
)

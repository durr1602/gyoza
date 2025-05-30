from snakemake.script import snakemake
from scripts.my_functions import load_codon_dic, annotate_mutants
import pandas as pd

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
        annot_df = pd.DataFrame(columns=[
            "nt_seq", "WT_seq",
            "WT", "aa_seq", "Nham_codons", "Nham_nt", "Nham_aa",
            "mutated_codon", "pos", "alt_codons", "alt_aa"
        ])

    # Ensure dataframe with indels is created even if empty
    if df_indels.empty:
        df_indels = pd.DataFrame(columns=df.columns)

    # Save outputs
    annot_df.to_csv(outpath, index=False)
    df_indels.to_csv(indel_outpath, index=False)

    return

get_annotated_mutants(snakemake.input[0],
                      snakemake.output.annot_rc,
                      snakemake.output.indels,
                      snakemake.config["codon"]["table"],
                      )

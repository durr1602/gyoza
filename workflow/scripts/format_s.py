from snakemake.script import snakemake
import pandas as pd
import seaborn as sns
import pickle


def get_heatmap_s_data(f, outpath, meta_out, pos_offset):
    """
    Takes dataframe of average functional impact scores as input.
    Each input dataframe should correspond to a group of samples with all selected time points.
    Formats df for each time point to be passed to the heatmap plotting module.
    Metadata are pickled.
    """

    AA_LIST = "*PGCQNTSEDKHRWYFMLIVA"
    AA_SORT = dict(zip(list(AA_LIST), list(range(0, len(AA_LIST)))))

    df = pd.read_csv(f)

    # Get min/max values
    vmax = max(1, int(df[[x for x in df.columns if "fitness_" in x]].max().max()) + 1)
    vmin = min(-1, int(df[[x for x in df.columns if "fitness_" in x]].min().min()))

    # Retrieve wild-type
    wtaa = df.loc[df.Nham_aa == 0, "aa_seq"].values[0]

    # Extract time point from outpath
    t = outpath.split("_format_s.csv")[0].split("_")[-1]

    # Duplicate WT for each position
    wt_fitness = df.loc[df.Nham_aa == 0, f"fitness_{t}"].values[0]
    wtdf = pd.DataFrame(
        [(int(i) + pos_offset, aa, wt_fitness) for i, aa in enumerate(wtaa)],
        columns=["aa_pos", "alt_aa", f"fitness_{t}"],
    )

    # Get single mutants
    singles = df[df.Nham_aa == 1][["aa_pos", "alt_aa", f"fitness_{t}"]].copy()
    singles["aa_pos"] = singles["aa_pos"].astype(int)

    # Assemble, pivot and sort dataframe
    filtered = pd.concat([wtdf, singles], ignore_index=True)
    wide = filtered.pivot(index="alt_aa", columns="aa_pos", values=f"fitness_{t}")
    wide.sort_index(key=lambda x: x.map(AA_SORT), inplace=True)

    # Export dataframe
    wide.to_csv(outpath)

    # Calculate WT coordinates
    wtcoord = [(i + 0.5, list(AA_LIST).index(v) + 0.5) for i, v in enumerate(wtaa)]

    # Set color map
    cmap = sns.color_palette(
        "blend:#009B9E,#42B7B9,#A7D3D4,#F1F1F1,#E4C1D9,#D691C1,#C75DAB",  # CARTOColors Tropic
        as_cmap=True,
    )
    cmap.set_bad(".5")

    # Assemble and save metadata
    meta = {
        "idx": "alt_aa",
        "fitness": f"fitness_{t}",
        "wt_coordinates": wtcoord,
        "color_map": cmap,
        "vmax": vmax,
        "vmin": vmin,
    }
    with open(meta_out, "wb") as m:
        pickle.dump(meta, m)

    return


get_heatmap_s_data(
    snakemake.input[0],
    snakemake.output.heatmap_df,
    snakemake.output.heatmap_meta,
    snakemake.params.position_offset,
)

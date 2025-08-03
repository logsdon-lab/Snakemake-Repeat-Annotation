import sys
import polars as pl


def main():
    # "/project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/"
    fai_dir = sys.argv[1]
    # "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/all.bed"
    bed = sys.argv[2]
    # Base pairs to add.
    try:
        bp_slop = int(sys.argv[3])
    except IndexError:
        bp_slop = 1_000_000

    df_fai = (
        pl.read_csv(f"{fai_dir}/*.fai", glob=True, has_header=False, separator="\t")
        .select("column_1", "column_2")
        .rename({"column_1": "chrom", "column_2": "chrom_length"})
    )
    df_bed = pl.read_csv(
        bed,
        separator="\t",
        new_columns=[
            "chrom",
            "chrom_st",
            "chrom_end",
            "div_time",
            "empty",
            "sample",
            "chrom_name",
            "clade_n",
            "chrom_og",
        ],
        has_header=False,
    )

    df_joined = (
        df_bed.join(df_fai, on="chrom", how="left")
        .with_columns(
            adj_chrom_st=pl.col("chrom_st") - bp_slop,
            adj_chrom_end=pl.col("chrom_end") + bp_slop,
        )
        .with_columns(
            adj_chrom_st_lt=pl.col("adj_chrom_st") < 0,
            adj_chrom_end_gt=pl.col("adj_chrom_end") > pl.col("chrom_length"),
        )
    )
    assert df_joined.filter(
        pl.col("adj_chrom_st_lt") | pl.col("adj_chrom_end_gt")
    ).is_empty(), "Slop exceeds contig boundaries."

    df_joined = df_joined.with_columns(
        chrom_st=pl.col("adj_chrom_st"),
        chrom_end=pl.col("adj_chrom_end"),
    ).drop(
        "adj_chrom_st",
        "adj_chrom_end",
        "adj_chrom_st_lt",
        "adj_chrom_end_gt",
        "chrom_length",
    )
    df_joined.write_csv(sys.stdout, separator="\t", include_header=False)


if __name__ == "__main__":
    main()

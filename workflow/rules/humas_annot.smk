HUMAS_ANNOT_OUTDIR = join(OUTPUT_DIR, "humas_annot")
HUMAS_ANNOT_LOGDIR = join(LOG_DIR, "humas_annot")
HUMAS_ANNOT_BMKDIR = join(BMK_DIR, "humas_annot")


if config["humas_annot"]["mode"] == "sd":
    humas_module_smk = "Snakemake-HumAS-SD/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-SD/workflow/envs/env.yaml"
else:
    humas_module_smk = "Snakemake-HumAS-HMMER/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-HMMER/workflow/envs/env.yaml"


module HumAS_Annot:
    snakefile:
        humas_module_smk
    config:
        {
            **config["humas_annot"],
            # if hmmer, otherwise no effect.
            "mode": "hor",
            "input_dir": join(SPLIT_MULTIFA_DIR),
            "output_dir": HUMAS_ANNOT_OUTDIR,
            "logs_dir": HUMAS_ANNOT_LOGDIR,
            "benchmarks_dir": HUMAS_ANNOT_BMKDIR,
        }


use rule * from HumAS_Annot as humas_*


rule convert_stv_bed_to_absolute_coords:
    input:
        bed=rules.humas_generate_stv.output,
    output:
        bed=join(HUMAS_ANNOT_OUTDIR, "{sm}_{fname}", "stv_row_abs.bed"),
    shell:
        """
        awk 'NR > 1 {{
            match($1, ":(.+)-", sts);
            $1 += sts[1]; $2 += sts[1];
            print
        }}' {input} > {output}
        """


# https://stackoverflow.com/a/63040288
def humas_annot_sm_outputs(wc):
    _ = checkpoints.split_multifasta.get(**wc).output
    wcs = glob_wildcards(
        join(SPLIT_MULTIFA_DIR, f"{wc.sm}_{{fname}}.fa"),
    )
    fnames = [f"{wc.sm}_{fname}" for fname in wcs.fname]

    return {
        "stv": expand(
            rules.convert_stv_bed_to_absolute_coords.output, zip, fname=fnames
        ),
    }


rule run_humas_annot:
    input:
        # # Force monomers to be generated.
        # (
        #     rules.humas_generate_monomers.output
        #     if config["humas_annot"]["mode"] == "sd"
        #     else []
        # ),
        expand(rules.split_multifasta.output, sm=SAMPLES),
        unpack(humas_annot_sm_outputs),
    output:
        touch(join(HUMAS_ANNOT_OUTDIR, "humas_annot_{sm}.done")),


rule humas_annot_all:
    input:
        expand(rules.run_humas_annot.output, sm=SAMPLES),
    default_target: True

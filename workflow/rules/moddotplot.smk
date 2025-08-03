MODDOTPLOT_OUTDIR = join(OUTPUT_DIR, "moddotplot")
MODDOTPLOT_LOGDIR = join(LOG_DIR, "moddotplot")
MODDOTPLOT_BMKDIR = join(BMK_DIR, "moddotplot")


rule moddotplot:
    input:
        fasta=join(SPLIT_MULTIFA_DIR, "{sm}_{fname}.fa"),
    output:
        plots=expand(
            join(
                MODDOTPLOT_OUTDIR,
                "{{sm}}",
                "{{fname}}",
                "{{fname}}_{otype}.{ext}",
            ),
            otype=["FULL", "HIST", "TRI"],
            ext=["png", "svg"],
        ),
        bed=join(MODDOTPLOT_OUTDIR, "{sm}", "{fname}", "{fname}.bed"),
    conda:
        "../envs/tools.yaml"
    params:
        window=config["moddotplot"]["window"],
        ident_thr=config["moddotplot"]["ident_thr"],
        outdir=lambda wc, output: os.path.dirname(output.bed),
    resources:
        mem=config["moddotplot"]["mem"],
    log:
        join(MODDOTPLOT_LOGDIR, "moddotplot_{sm}_{fname}.log"),
    benchmark:
        join(MODDOTPLOT_BMKDIR, "moddotplot_{sm}_{fname}.tsv")
    # singularity:
    #     "/project/logsdon_shared/tools/moddotplot.sif"
    shell:
        """
        moddotplot static -f {input.fasta} -w {params.window} -o {params.outdir} -id {params.ident_thr} &> {log}
        """


# Gather all RM output
def moddotplot_output(wc):
    _ = checkpoints.split_multifasta.get(**wc).output
    fnames = glob_wildcards(join(SPLIT_MULTIFA_DIR, f"{wc.sm}_{{fname}}.fa")).fname
    return expand(rules.moddotplot.output, sm=wc.sm, fname=fnames)


rule run_moddotplot:
    input:
        moddotplot_output,
    output:
        touch(join(MODDOTPLOT_OUTDIR, "{sm}.done")),


rule moddotplot_all:
    input:
        expand(rules.run_moddotplot.output, sm=SAMPLES),

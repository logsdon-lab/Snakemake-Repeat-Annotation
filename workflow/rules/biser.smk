BISER_OUTDIR = join(OUTPUT_DIR, "biser")
BISER_LOGDIR = join(LOG_DIR, "biser")
BISER_BMKDIR = join(BMK_DIR, "biser")


rule run_biser:
    input:
        fa=rules.repeatmasker_softmasked_assembly.output.seq,
    output:
        join(BISER_OUTDIR, "{sm}.bedpe"),
    log:
        join(BISER_LOGDIR, "biser_{sm}.log"),
    benchmark:
        join(BISER_BMKDIR, "biser_{sm}.tsv")
    conda:
        "../envs/tools.yaml"
    threads: config["biser"]["threads"]
    resources:
        mem=config["biser"]["mem"],
        gc_heap=config["biser"]["mem"].strip("B"),
    shell:
        """
        # https://github.com/0xTCG/biser/issues/32
        biser -o {output} -t {threads} --keep-contigs --gc-heap {resources.gc_heap} {input.fa}
        """


rule biser_all:
    input:
        expand(rules.run_biser.output, sm=SAMPLES),

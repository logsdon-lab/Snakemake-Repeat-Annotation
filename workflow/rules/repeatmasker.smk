RM_OUTDIR = join(OUTPUT_DIR, "repeatmasker")
RM_LOGDIR = join(LOG_DIR, "repeatmasker")
RM_BMKDIR = join(BMK_DIR, "repeatmasker")


rule setup_repeatmasker:
    output:
        chkpt=touch(join(RM_OUTDIR, "rm_setup.done")),
        seq=join(RM_OUTDIR, "rm_setup.fa"),
        rm_dir=directory(join(RM_OUTDIR, "rm_setup")),
    params:
        species=config["repeatmasker"]["species"],
        engine=config["repeatmasker"]["engine"],
    threads: 1
    log:
        join(RM_LOGDIR, "setup_repeatmasker.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        echo ">rm_setup" > {output.seq}
        echo "NNNNNNNNNNNNNNNNNNNNN" >> {output.seq}
        RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {output.rm_dir} \
            -pa {threads} \
            {output.seq} &> {log}
        """


rule rename_for_repeatmasker:
    input:
        fa=join(SPLIT_MULTIFA_DIR, "{fname}.fa"),
    output:
        original_fa_idx=temp(
            join(RM_OUTDIR, "fa", "{sm}_renamed", "{fname}_original.fa.fai"),
        ),
        renamed_fa=temp(
            join(
                RM_OUTDIR,
                "seq",
                "{sm}_renamed",
                "{fname}.fa",
            )
        ),
        renamed_fa_idx=temp(
            join(
                RM_OUTDIR,
                "seq",
                "{sm}_renamed",
                "{fname}.fa.fai",
            )
        ),
    params:
        prefix="seq",
    conda:
        "../envs/tools.yaml"
    log:
        join(RM_LOGDIR, "rename_for_repeatmasker_{sm}_{fname}.log"),
    shell:
        """
        samtools faidx {input.fa} -o {output.original_fa_idx} 2> {log}
        seqtk rename {input.fa} {params.prefix} > {output.renamed_fa} 2>> {log}
        if [ -s {output.renamed_fa} ]; then
            samtools faidx {output.renamed_fa} 2>> {log}
        else
            touch {output.renamed_fa_idx}
        fi
        """


rule run_repeatmasker:
    input:
        setup=rules.setup_repeatmasker.output,
        seq=rules.rename_for_repeatmasker.output.renamed_fa,
    output:
        temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}_renamed",
                "{fname}.fa.out",
            )
        ),
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: dirname(str(output)),
        species=config["repeatmasker"]["species"],
        engine=config["repeatmasker"]["engine"],
    resources:
        mem=config["repeatmasker"]["mem"],
    conda:
        "../envs/tools.yaml"
    log:
        join(RM_LOGDIR, "repeatmasker_{sm}_{fname}.log"),
    benchmark:
        join(RM_BMKDIR, "repeatmasker_{sm}_{fname}.tsv")
    shell:
        """
        RepeatMasker \
        -engine {params.engine} \
        -species {params.species} \
        -dir {params.output_dir} \
        -pa {threads} \
        {input.seq} &> {log}
        """


# Rename repeatmasker output to match the original sequence names.
rule reformat_repeatmasker_output:
    input:
        rm_out=rules.run_repeatmasker.output,
        original_fai=rules.rename_for_repeatmasker.output.original_fa_idx,
        renamed_fai=rules.rename_for_repeatmasker.output.renamed_fa_idx,
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "{sm}",
            "{fname}.fa.out",
        ),
    params:
        script=workflow.source_path("../scripts/rename_rm.py"),
    log:
        join(RM_LOGDIR, "reformat_repeatmasker_output_{sm}_{fname}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        python {params.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} > {output} 2> {log}
        """


# Gather all RM output
def refmt_rm_output(wc):
    _ = checkpoints.split_multifasta.get(**wc).output
    fa_glob_pattern = join(SPLIT_MULTIFA_DIR, "{fname}.fa")
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    return expand(rules.reformat_repeatmasker_output.output, sm=wc.sm, fname=fnames)


rule repeatmasker_output:
    input:
        rm_out=refmt_rm_output,
        # Force snakemake to not evaluate chkpt function until all dirs created.
        rm_fa_chkpt=expand(rules.split_multifasta.output, sm=SAMPLES),
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "all",
            "{sm}_cens.fa.out",
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_out} > {output}
        """


rule repeatmasker_all:
    input:
        expand(rules.repeatmasker_output.output, sm=SAMPLES),
    default_target: True

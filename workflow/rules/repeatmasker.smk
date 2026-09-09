RM_OUTDIR = join(OUTPUT_DIR, "repeatmasker")
RM_LOGDIR = join(LOG_DIR, "repeatmasker")
RM_BMKDIR = join(BMK_DIR, "repeatmasker")


rule setup_repeatmasker:
    output:
        chkpt=touch(join(RM_OUTDIR, "rm_setup.done")),
        seq=join(RM_OUTDIR, "rm_setup.fa"),
        rm_dir=directory(join(RM_OUTDIR, "rm_setup")),
    log:
        join(RM_LOGDIR, "setup_repeatmasker.log"),
    conda:
        "../envs/tools.yaml"
    threads: 1
    params:
        species=config["repeatmasker"]["species"],
        engine=config["repeatmasker"]["engine"],
    shell:
        """
        echo ">rm_setup" >{output.seq}
        echo "NNNNNNNNNNNNNNNNNNNNN" >>{output.seq}
        RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {output.rm_dir} \
            -pa {threads} \
            {output.seq} &>{log}
        """


rule rename_for_repeatmasker:
    input:
        fa=join(SPLIT_MULTIFA_DIR, "{sm}_{fname}.fa"),
    output:
        original_fa_idx=temp(
            join(RM_OUTDIR, "fa", "{sm}_original", "{fname}.fa.fai"),
        ),
        renamed_fa=temp(
            join(
                RM_OUTDIR,
                "fa",
                "{sm}_renamed",
                "{fname}.fa",
            )
        ),
        renamed_fa_idx=temp(
            join(
                RM_OUTDIR,
                "fa",
                "{sm}_renamed",
                "{fname}.fa.fai",
            )
        ),
    log:
        join(RM_LOGDIR, "rename_for_repeatmasker_{sm}_{fname}.log"),
    conda:
        "../envs/tools.yaml"
    params:
        prefix="seq",
    shell:
        """
        samtools faidx {input.fa} -o {output.original_fa_idx} 2>{log}
        seqtk rename {input.fa} {params.prefix} >{output.renamed_fa} 2>>{log}
        if [ -s {output.renamed_fa} ]; then
            samtools faidx {output.renamed_fa} 2>>{log}
        else
            touch {output.renamed_fa_idx}
        fi
        """


rule run_repeatmasker:
    input:
        setup=rules.setup_repeatmasker.output,
        seq=rules.rename_for_repeatmasker.output.renamed_fa,
    output:
        rm_seq=temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}_renamed",
                "{fname}.fa.masked",
            )
        ),
        rm_out=temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}_renamed",
                "{fname}.fa.out",
            )
        ),
    log:
        join(RM_LOGDIR, "repeatmasker_{sm}_{fname}.log"),
    benchmark:
        join(RM_BMKDIR, "repeatmasker_{sm}_{fname}.tsv")
    conda:
        "../envs/tools.yaml"
    threads: config["repeatmasker"]["threads"]
    resources:
        mem=config["repeatmasker"]["mem"],
    params:
        output_dir=lambda wc, output: dirname(str(output.rm_out)),
        species=config["repeatmasker"]["species"],
        engine=config["repeatmasker"]["engine"],
    shell:
        """
        {{ RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {params.output_dir} \
            -xsmall \
            -pa {threads} \
            {input.seq} || true ;}} &>{log}
        touch {output}
        """


# Rename repeatmasker output to match the original sequence names.
# Also convert to original coordinate system if start in name
rule reformat_repeatmasker_output:
    input:
        rm_seq=rules.run_repeatmasker.output.rm_seq,
        rm_out=rules.run_repeatmasker.output.rm_out,
        original_fai=rules.rename_for_repeatmasker.output.original_fa_idx,
        renamed_fai=rules.rename_for_repeatmasker.output.renamed_fa_idx,
    output:
        rm_seq=temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}",
                "{fname}.fa.masked",
            )
        ),
        rm_out=temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}",
                "{fname}.fa.out",
            )
        ),
    log:
        join(RM_LOGDIR, "reformat_repeatmasker_output_{sm}_{fname}.log"),
    conda:
        "../envs/tools.yaml"
    params:
        script=workflow.source_path("../scripts/rename_rm.py"),
    shell:
        """
        python {params.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} >{output.rm_out} 2>{log}
        awk '{{ if ($1 ~ ">") {{ $1=">{wildcards.fname}" }} print }}' {input.rm_seq} >{output.rm_seq}
        """


# Gather all RM output
def refmt_rm_output(wc, typ: str):
    _ = checkpoints.split_multifasta.get(**wc).output
    fnames = glob_wildcards(join(SPLIT_MULTIFA_DIR, f"{wc.sm}_{{fname}}.fa")).fname

    return expand(
        getattr(rules.reformat_repeatmasker_output.output, typ), sm=wc.sm, fname=fnames
    )


rule repeatmasker_output:
    input:
        rm_out=lambda wc: refmt_rm_output(wc, "rm_out"),
        # Force snakemake to not evaluate chkpt function until all dirs created.
        rm_fa_chkpt=expand(rules.split_multifasta.output, sm=SAMPLES),
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "all",
            "{sm}.fa.out",
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_out} >{output}
        """


rule repeatmasker_softmasked_assembly:
    input:
        rm_seq=lambda wc: refmt_rm_output(wc, "rm_seq"),
        # Force snakemake to not evaluate chkpt function until all dirs created.
        rm_fa_chkpt=expand(rules.split_multifasta.output, sm=SAMPLES),
    output:
        seq=join(
            RM_OUTDIR,
            "fa",
            "all",
            "{sm}.fa.gz",
        ),
        idx=join(
            RM_OUTDIR,
            "fa",
            "all",
            "{sm}.fa.gz.fai",
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        cat {input.rm_seq} | bgzip >{output.seq}
        samtools faidx {output.seq}
        """


rule repeatmasker_all:
    default_target: True
    input:
        expand(rules.repeatmasker_output.output, sm=SAMPLES),

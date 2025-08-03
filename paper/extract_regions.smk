import os


ASM_DIR = "/project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/"
with open(
    "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/all.bed"
) as fh:
    SAMPLES = []
    for line in fh:
        elems = line.strip().split("\t")
        SAMPLES.append(elems[5])


wildcard_constraints:
    sm="|".join(SAMPLES),


rule slop_bed:
    input:
        bed="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/all.bed",
        fa_dir=ASM_DIR,
    output:
        bed="results/bed/all_1mbp_slop.bed",
    params:
        script=workflow.source_path("slop_bed.py"),
        bp_slop=1_000_000,
    shell:
        """
        python {params.script} {input.fa_dir} {input.bed} {params.bp_slop} > {output}
        """


rule extract_regions:
    input:
        bed=rules.slop_bed.output,
        asm=os.path.join(ASM_DIR, "{sm}-asm-renamed-reort.fa"),
    output:
        "results/split_fa/{sm}.fa.gz",
    shell:
        """
        seqtk subseq {input.asm} <(cut -f 1-3 {input.bed}) | bgzip > {output}
        samtools faidx {output}
        """


rule make_config:
    input:
        fa=expand(rules.extract_regions.output, sm=SAMPLES),
        cfg="/project/logsdon_shared/projects/HGSVC3/Snakemake-Repeat-Annotation/config/config.yaml",
    output:
        "results/cfg/annot.yaml",
    params:
        script=workflow.source_path("make_cfg.py"),
    shell:
        """
        python {params.script} -c {input.cfg} -i {input.fa} > {output}
        """


rule all:
    input:
        rules.make_config.output,
    default_target: True

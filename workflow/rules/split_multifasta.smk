SPLIT_MULTIFA_DIR = join(OUTPUT_DIR, "fa")


checkpoint split_multifasta:
    input:
        fa=lambda wc: config["samples"][wc.sm]["fa"],
        bed=lambda wc: config["samples"][wc.sm].get("bed", []),
    output:
        join(
            SPLIT_MULTIFA_DIR,
            "{sm}.fofn",
        ),
    log:
        join(LOG_DIR, "split_multifasta_{sm}.log"),
    params:
        output_dir=SPLIT_MULTIFA_DIR,
        extract_region=lambda wc, input: (
            f"| seqtk subseq - {input.bed}" if input.bed else ""
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.output_dir}
        awk '{{
            if (substr($0, 1, 1)==">") {{
                filename=("{params.output_dir}/{wildcards.sm}_" substr($0,2) ".fa")
                print filename
            }}
            print $0 > filename
        }}' <(zcat -f {input.fa} {params.extract_region}) > {output} 2> {log}
        """

################################################################################
# Project: Streptococcus mutans gene/presence analysis
# Part:
#
#	NOTE: THIS IS COMMAND LINE FORMAT: 
#	snakemake --configfile ./probes_pangenome_mapq0/gene_content_analysis-config.json --cores 4 -s gene_content_analysis.smk 
#
#
# Megan Michel, 01/30/2020
################################################################################

import pandas as pd
workdir: config['workdir']

SAMPLES, = glob_wildcards(config['BAM_DIRECTORY'] + "{sample}.bam")
print (SAMPLES)

rule all:
    input:
        expand("{sample}.outfile.bed", sample=SAMPLES)

rule get_genome_cov:
     output:
         temp("{sample}.histogram.bed")
     params:
         config['BAM_DIRECTORY'] + "{sample}.bam"
     shell:
         """
         bedtools genomecov -ibam {params} -bg > {output}
         """

rule a:
    input:
        sample_histo = "{sample}.histogram.bed",
        target_genes = config['TARGET_GENE_BED']
    output:
        vir1 = temp("{sample}.virulence_only1.bed"),
        vir0 = temp("{sample}.virulence_only0.bed")
    shell:
        """
        bedtools coverage -a {input.target_genes} -b {input.sample_histo} -hist \
        | awk '$4==0' | awk '$7==1.0000000'  > {output.vir0}
        bedtools coverage -a {input.target_genes} -b {input.sample_histo} -hist \
        | awk '$4==1' > {output.vir1}
        """

rule b:
    input:
        "{sample}.virulence_only0.bed"
    output:
        temp("{sample}.virulence_only0_with0.bed")
    shell:
        """
        sed 's/\t1.0000000/\t0.0000000/g' {input} > {output}
        """

rule c:
    input:
        only0 = "{sample}.virulence_only0_with0.bed",
        only1 = "{sample}.virulence_only1.bed"
    output:
        temp("{sample}.virulence.bed")
    shell:
        """
        cat {input.only1} {input.only0} > {output}
        """

rule e:
    input:
        "{sample}.virulence.bed"
    output:
        temp("{sample}.virulence.sorted.bed")
    shell:
        """
        bedtools sort -i {input} > {output}
        """

rule f:
    input:
        "{sample}.virulence.sorted.bed"
    output:
        temp("{sample}.virulence.sorted.filtered.bed")
    params:
        name = "{sample}"
    shell:
        """
        cat {input} | cut -f 1,7 |awk '{{$1="$NAME"; print ;}}' \
        | sed 's/\$NAME/{params.name}/g' > {output}
        """

rule g:
    input:
        target_gene_names = config['TARGET_GENE_NAMES'],
        virulence_sorted_filtered = "{sample}.virulence.sorted.filtered.bed"
    output:
        "{sample}.outfile.bed"
    shell:
        """
        paste <(awk '{{print $1 "\t" $2}}' {input.virulence_sorted_filtered} )\
         <(awk -F "\t" '{{print $1}}' {input.target_gene_names}) > {output}
        """



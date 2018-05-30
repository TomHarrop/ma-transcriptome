#!/usr/bin/env python3

import pathlib2


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))


def separate_input_reads(wildcards):
    if wildcards.sample == 's1':
        return {
            'r1': 'data/reads/C7NA2ANXX-1865-07-1-1_L006_R1.fastq',
            'r2': 'data/reads/C7NA2ANXX-1865-07-1-1_L006_R2.fastq'}
    elif wildcards.sample == 's2':
        return {
            'r1': 'data/reads/C7NA2ANXX-1865-08-1-1_L006_R1.fastq',
            'r2': 'data/reads/C7NA2ANXX-1865-08-1-1_L006_R2.fastq'}
    else:
        raise ValueError('Sample {} not found'.format(wildcards.sample))

###########
# GLOBALS #
###########

read_dir = 'data/reads'
bbduk_adaptors = '/adapters.fa'

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
trinity_container = 'shub://TomHarrop/singularity-containers:trinity_2.6.6'


#########
# RULES #
#########

rule target:
    input:
        expand('output/busco/run_{filter}/full_table_all.tsv',
               filter=['expression', 'length'])

wildcard_constraints:
    sample = 's\d'

rule busco:
    input:
        fasta = 'output/filtered_isoforms/isoforms_by_{filter}.fasta',
        lineage = 'data/lineages/hymenoptera_odb9'
    output:
        'output/busco/run_{filter}/full_table_all.tsv'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                          'busco_{filter}.log'))
    benchmark:
        'output/benchmark/busco_{filter}.tsv'
    params:
        wd = 'output/busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage)
    threads:
        20
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--in {params.fasta} '
        '--out {wildcards.filter} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species honeybee1 '
        '--mode transcriptome '
        '&> {log}'

rule filter_trinity_isoforms:
    input:
        fasta = 'output/trinity/Trinity.fasta',
        isoforms = 'output/filtered_isoforms/isoform_by_{filter}.txt'
    output:
        fasta = 'output/filtered_isoforms/isoforms_by_{filter}.fasta'
    log:
        'output/logs/filter_isoforms_by_{filter}.log'
    benchmark:
        'output/benchmark/filter_isoforms_by_{filter}.tsv'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input.fasta} '
        'include=t '
        'names={input.isoforms} '
        'out={output.fasta} '
        '&> {log}'

rule select_trinity_isoforms:
    input:
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        expression = 'output/filtered_isoforms/isoform_by_expression.txt',
        length = 'output/filtered_isoforms/isoform_by_length.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/select_trinity_isoforms.log'
    benchmark:
        'output/benchmark/select_trinity_isoforms.tsv'
    script:
        'src/select_trinity_isoforms.R'

rule trinity_abundance_to_matrix:
    input:
        gt_map = 'output/trinity/Trinity.fasta.gene_trans_map',
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        'output/trinity_abundance/RSEM.isoform.counts.matrix',
        'output/trinity_abundance/RSEM.isoforms.results'
    params:
        prefix = 'output/trinity_abundance/RSEM'
    singularity:
        trinity_container
    log:
        'output/logs/abundance_estimates_to_matrix.log'
    benchmark:
        'output/benchmark/abundance_estimates_to_matrix.tsv'
    shell:
        'abundance_estimates_to_matrix.pl '
        '--est_method RSEM '
        '--cross_sample_norm none '
        '--out_prefix {params.prefix} '
        '--gene_trans_map {input.gt_map} '
        '{input.abundance} '
        '&> {log}'

rule trinity_abundance:
    input:
        transcripts = 'output/trinity/Trinity.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz',
                      sample=['s1', 's2']),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz',
                       sample=['s1', 's2'])
    output:
        'output/trinity_abundance/RSEM.isoforms.results'
    params:
        outdir = 'output/trinity_abundance',
        left = lambda wilcdards, input: ','.join(input.left),
        right = lambda wilcdards, input: ','.join(input.right)
    singularity:
        trinity_container
    threads:
        20
    log:
        'output/logs/align_and_estimate_abundance.log'
    benchmark:
        'output/benchmark/align_and_estimate_abundance.tsv'
    shell:
        'align_and_estimate_abundance.pl '
        '--transcripts {input.transcripts} '
        '--seqType fq '
        '--est_method RSEM '
        '--output_dir {params.outdir} '
        '--aln_method bowtie '
        '--prep_reference '
        '--SS_lib_type RF '
        '--thread_count 20 '
        '--trinity_mode '
        '--left {params.left} '
        '--right {params.right} '
        '&> {log}'


rule assemble_transcriptome:
    input:
        left = expand('output/bbmerge/{sample}_r1_all.fq.gz',
                      sample=['s1', 's2']),
        right = expand('output/bbmerge/{sample}_r2_unmerged.fq.gz',
                       sample=['s1', 's2'])
    output:
        'output/trinity/Trinity.fasta',
        'output/trinity/Trinity.fasta.gene_trans_map'
    params:
        outdir = 'output/trinity',
        left = lambda wilcdards, input: ','.join(input.left),
        right = lambda wilcdards, input: ','.join(input.right) 
    singularity:
        trinity_container
    threads:
        20
    log:
        'output/logs/trinity.log'
    benchmark:
        'output/benchmark/trinity.tsv'
    shell:
        'Trinity '
        '--SS_lib_type RF '
        '--max_memory 500G '
        '--CPU {threads} '
        '--output {params.outdir} '
        '--left {params.left} '
        '--right {params.right} '
        '--seqType fq '
        '&> {log}'


rule add_merged_to_r1_file:
    input:
        r1 = 'output/bbmerge/{sample}_r1_unmerged.fq.gz',
        merged = 'output/bbmerge/{sample}_merged.fq.gz'
    output:
        joined = 'output/bbmerge/{sample}_r1_all.fq.gz'
    shell:
        'cat {input.r1} {input.merged} > {output.joined}'

rule bbmerge:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        merged = 'output/bbmerge/{sample}_merged.fq.gz',
        ur1 = 'output/bbmerge/{sample}_r1_unmerged.fq.gz',
        ur2 = 'output/bbmerge/{sample}_r2_unmerged.fq.gz',
        ihist = 'output/bbmerge/{sample}_ihist.txt'
    params:
        ref = bbduk_adaptors
    log:
        'output/logs/bbmerge_{sample}.log'
    benchmark:
        'output/benchmark/bbmerge_{sample}.tsv'
    singularity:
        bbduk_container
    threads:
        20
    shell:
        'bbmerge.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.merged} '
        'outu1={output.ur1} '
        'outu2={output.ur2} '
        'ihist={output.ihist} '
        'verystrict=t '
        'adapters={params.ref} '
        '&> {log}'

rule bbduk_trim:
    input:
        unpack(separate_input_reads)
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        ref = bbduk_adaptors
    log:
        'output/logs/bbduk_{sample}.log'
    benchmark:
        'output/benchmark/bbduk_{sample}.tsv'
    singularity:
        bbduk_container
    threads:
        20
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.ref} '
        'ktrim=r '
        'k=23 '
        'mink=11 '
        'hdist=1 '
        'tpe tbo '
        'qtrim=r '
        'trimq=15 '
        '&> {log}'


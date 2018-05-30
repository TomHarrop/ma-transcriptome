#!/usr/bin/env python3


#############
# FUNCTIONS #
#############

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
bbduk_container = ('shub:// TomHarrop/singularity-containers:'
                   'bbmap_38.00')
trinity_container = ('shub:// TomHarrop/singularity-containers:'
                     'trinity_2.6.6')


#########
# RULES #
#########

rule target:
    input:
        'output/trinity_abundance/RSEM.isoform.counts.matrix'

wildcard_constraints:
    sample = 's\d'

rule trinity_abundance_to_matrix:
    input:
        gt_map = 'output/trinity/Trinity.fasta.gene_trans_map',
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        'output/trinity_abundance/RSEM.isoform.counts.matrix'
    params:
        prefix = 'output/trinity_abundance/RSEM'
    singularity:
        trinity_container
    log:
        'output/logs/abundance_estimates_to_matrix.log'
    shell:
        'abundance_estimates_to_matrix.pl '
        '--est_method RSEM '
        '--cross_sample_norm none '
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
    singularity:
        bbduk_container
    shell:
        'bbmerge.sh '
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
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.ref} '
        'ktrim=5 '
        'k=23 '
        'mink=11 '
        'hdist=1 '
        'tpe tbo '
        'qtrim=r '
        'trimq=15 '
        '&> {log}'


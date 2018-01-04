#!/usr/bin/env snakemake

"""
Snakemake QTL pipeline
----------------------

Snakemake pipeline for molecular trait QTL mapping
"""

__version__ = '0.1.4.dev0'

rule all:
    input:
        # QTLtools nominal associations
        'qtltools_nominal.tsv.gz',
        
        # QTLtools permutation associations
        'qtltools_permute.tsv.gz',
        'plots/qtltools_permute-check_beta_approx.png', 
        
        # QTLtools conditional associations
        'qtltools_conditional.tsv.gz', 
        
        # veqtl-mapper variance qtls
        'vqtlmapper_permute-top_qtl.tsv.gz'


# QTLtools: file prep ##########################################################
rule qtltools_prep__vcf:
    """
    QTLtools:   generate tabix index of vcf input files 
    """
    input:
        'data/genotypes.vcf.gz'
    output:
        'data/genotypes.vcf.gz.tbi'
    shell:
        'tabix {input}'


rule qtltools_prep__bed:
    """
    QTLtools:   generate tabix index of bed input files 
    """
    input:
        'data/moltraits.bed.gz'
    output:
        'data/moltraits.bed.gz.tbi'
    shell:
        'tabix {input}'
################################################################################


# QTLtools: nominal pass #######################################################
rule qtltools_cis__nominal:
    """
    QTLtools:   nominal cis mapping (association with all variants)
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits.bed.gz',
        pheno_tbi='data/moltraits.bed.gz.tbi'
    output:
        'swarm_qtltools_nominal/qtltools_nominal-{j_cur}_{j_total}.txt'
    params:
        options=config['OPTIONS_QTLTOOLS']
    run:
        shell(
            'QTLtools cis --vcf {input.geno} --bed {input.pheno} '
                '--nominal 1 ' # nominal p-value cutoff to threshold output  
                '--normal ' # force the input phenotype to N(0, 1) 
                '--seed 15112011 ' # seed for random number generator
                '--window 1000000 ' # window size. default = 1Mb
                '--chunk {wildcards.j_cur} {wildcards.j_total} ' # for swarming
                '{params.options} ' # options from config
                '--out {output}'
        )

rule qtltools_cis__nominal_concat:
    """
    QTLtools: swarm cis mapping & combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    - removes the temp input data
    """
    input:
        expand(
            'swarm_qtltools_nominal/qtltools_nominal-{j_cur}_{j_total}.txt', 
            j_cur=range(1, config['NJOBS_NOMINAL']+1), 
            j_total=config['NJOBS_NOMINAL']
        )
    output:
        'qtltools_nominal.tsv.gz'
    run:
        # add the header to the output file
        shell(
            'echo "pheno_id '
                'pheno_chr pheno_start pheno_end pheno_strand '
                'n_cis_var pheno_var_dist '
                'var_id var_chr var_start var_end '
                'p_nominal beta top_cis_var" | '
            'sed s/" "/"\t"/g | gzip -c > {output}'
        )
            
        # cat swarm files to make the output file
        shell('cat {input} | sed s/" "/"\t"/g | gzip -c >> {output}')
        
        # remove the temp files and swarm directory
        shell('rm {input}')
        shell('rm -r swarm_qtltools_nominal')
################################################################################


# QTLtools: permute pass #######################################################
rule qtltools_cis__permute:
    """
    QTLtools:   cis mapping where genotypes are permuted and top cis variant
                returned
                
    Notes:
    - only cases where the first qtl is < threshold are in the output
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits.bed.gz',
        pheno_tbi='data/moltraits.bed.gz.tbi'
    output:
        'swarm_qtltools_permute/qtltools_permute-{j_cur}_{j_total}.txt'
    params:
        options=config['OPTIONS_QTLTOOLS'],
	n_perm=config['NPERM']
    run:
        shell(
            'QTLtools cis --vcf {input.geno} --bed {input.pheno} '
                '--permute {params.n_perm} ' # qtltools not have adaptive perm
                '--normal ' # force the input phenotype to N(0, 1) 
                '--seed 15112011 ' # seed for random number generator
                '--window 1000000 ' # window size. default = 1Mb
                '--chunk {wildcards.j_cur} {wildcards.j_total} ' # for swarming
                '{params.options} ' # options from config
                '--out {output}'
        )

rule qtltools_cis__permute_concat:
    """
    QTLtools:   swarm cis mapping & combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    - removes the temp input data
    """
    input:
        expand(
            'swarm_qtltools_permute/qtltools_permute-{j_cur}_{j_total}.txt',
            j_cur=range(1, config['NJOBS_PERMUTE']+1), 
            j_total=config['NJOBS_PERMUTE']
        )
    output:
        'qtltools_permute.tsv.gz'
    run:
        # add the header to the output file
        shell(
            'echo "pheno_id '
                'pheno_chr pheno_start pheno_end pheno_strand '
                'n_cis_var pheno_var_dist '
                'var_id var_chr var_start var_end '
                'p_degree_freedom '
                'dummy '
                'beta_dist_par1 beta_dist_par2 '
                'p_nominal beta '
                'p_empirical p_adjust_beta_dist" | '
            'sed s/" "/"\t"/g | gzip -c > {output}'
        )
            
        # cat swarm files to make the output file
        shell('cat {input} | sed s/" "/"\t"/g | gzip -c >> {output}')
        
        # remove the temp swarm directory
        shell('rm {input}')
        shell('rm -r swarm_qtltools_permute')

rule qtltools_cis__permute_fdr:
    """
    QTLtools:   run the FDR script for the permuted data
    """
    input:
        'qtltools_permute.tsv.gz'
    output:
        'qtltools_permute-significant.tsv.gz',
        'qtltools_permute-thresholds.txt.gz'
    params:
        shell_script=srcdir('scripts/qtltools-runFDR_cis.R')
    shell:
        'Rscript {params.shell_script} {input} 0.05 `basename {input} .tsv.gz`'

rule qtltools_cis__plot__check_beta_approx:
    input:
        'qtltools_permute.tsv.gz'
    output:
        'plots/qtltools_permute-check_beta_approx.png',
        'plots/qtltools_permute-check_beta_approx-neglog10.png'
    params:
        shell_script=srcdir('scripts/qtltools-check_beta_approx.R')
    shell:
        'Rscript {params.shell_script} {input} '
        'plots/`basename {input} .tsv.gz`-check_beta_approx'
################################################################################


# QTLtools: conditional pass ###################################################
rule qtltools_cis__conditional_unzip_threshold:
    """
    QTLtools requires uncompressed thresholds file.
    """
    input:
        'qtltools_permute-thresholds.txt.gz'
    output:
        'qtltools_permute-thresholds.txt'
    run:
        shell('gunzip -c {input} > {output}')
        
rule qtltools_cis__conditional:
    """
    QTLtools:   conditional cis mapping using thresholds from permuted data
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits.bed.gz',
        pheno_tbi='data/moltraits.bed.gz.tbi',
        pheno_thresh='qtltools_permute-thresholds.txt'
    output:
        'swarm_qtltools_conditional/qtltools_conditional-{j_cur}_{j_total}.txt'
    params:
        options=config['OPTIONS_QTLTOOLS']
    run:
        shell(
            'QTLtools cis --vcf {input.geno} --bed {input.pheno} '
                '--mapping {input.pheno_thresh} ' # thresholds
                '--normal ' # force the input phenotype to N(0, 1) 
                '--seed 15112011 ' # seed for random number generator
                '--window 1000000 ' # window size. default = 1Mb
                '--chunk {wildcards.j_cur} {wildcards.j_total} ' # for swarming
                '{params.options} ' # options from config 
                '--out {output}'
        )

rule qtltools_cis__conditional_concat:
    """
    QTLtools:   swarm cis mapping & combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    - removes the temp input data
    """
    input:
        expand(
            'swarm_qtltools_conditional/qtltools_conditional-{j_cur}_{j_total}.txt',
            j_cur=range(1, config['NJOBS_CONDITIONAL']+1), 
            j_total=config['NJOBS_CONDITIONAL']
        )
    output:
        'qtltools_conditional.tsv.gz'
    run:
        # add the header to the output file
        shell(
            'echo "pheno_id '
                'pheno_chr pheno_start pheno_end pheno_strand '
                'n_cis_var pheno_var_dist '
                'var_id var_chr var_start var_end '
                'var_rank ' # if variant = 1st best (rank=0), 2nd best (rank=1)
                'p_nominal_forward beta_forward '
                'var_top_rank_forward var_p_below_threshold_forward '
                'p_nominal_backward beta_backward '
                'var_top_rank_backward var_p_below_threshold_backward" | '
            'sed s/" "/"\t"/g | gzip -c > qtltools_conditional.tsv.gz'
        )
        
        # cat swarm files to make the output file
        shell(
            'cat {input} | sed s/" "/"\t"/g | '
            'gzip -c >> qtltools_conditional.tsv.gz'
        )
        
        # compress the thresholds file
        shell('gzip -f qtltools_permute-thresholds.txt')
        
        # remove the temp swarm directory
        shell('rm {input}')
        shell('rm -r swarm_qtltools_conditional')
################################################################################


# vqtl-mapper: file prep ######################################################
rule vqtlmapper_prep__bed_gz:
    """
    vqtl-mapper:    change bed input file to be consistent with fastqtl and 
                    vqtl-mapper
    """
    input:
        'data/moltraits.bed.gz'
    output:
        'data/moltraits-fastqtl.bed.gz'
    shell:
        'gunzip -c {input} | cut -f 1-4,7- | bgzip -c > '   
            'data/moltraits-fastqtl.bed.gz; '

rule vqtlmapper_prep__bed:
    """
    vqtl-mapper:    change bed input file to be consistent with fastqtl and 
                    vqtl-mapper
    """
    input:
        'data/moltraits-fastqtl.bed.gz'
    output:
        'data/moltraits-fastqtl.bed'
    shell:
        'gunzip -c {input} > data/moltraits-fastqtl.bed'
        
rule vqtlmapper_prep__bed_tabix:
    """
    vqtl-mapper:    change bed input file to be consistent with fastqtl and 
                    vqtl-mapper
    """
    input:
        'data/moltraits-fastqtl.bed.gz'
    output:
        'data/moltraits-fastqtl.bed.gz.tbi'
    shell:
        'tabix {input}'
################################################################################

# vqtl-mapper: permutation pass ###############################################
rule vqtlmapper__permute:
    """
    vqtl-mapper:   map variance qtls
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits-fastqtl.bed'
    output:
        'swarm_vqtl_permute/vqtl_permute-{j_cur}_{j_total}.txt'
    params:
        perm=config['NPERM'],
        seed='15112011',
        window='1000000',
        other_settings='--normal'
    script:
        'wrappers/veqtl-mapper/wrapper.py'

rule vqtlmapper__permute_concat:
    """
    vqtl-mapper:   swarm cis mapping & combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    - removes the temp input data
    - corrected for multiple testing across cis variants per feature
    """
    input:
        expand('swarm_vqtl_permute/vqtl_permute-{j_cur}_{j_total}.txt',
            j_cur=range(1, config['NJOBS_PERMUTE']+1), 
            j_total=config['NJOBS_PERMUTE'])
    output:
        'vqtlmapper_permute.tsv.gz'
    run:            
        # cat swarm files to make the output file
        shell(
            'awk "FNR>1||NR==1" {input} | sed s/"GENE"/"pheno_id"/ | '
            'gzip -c > {output}'
        )
        
        # remove the temp swarm directory
        shell('rm {input}')
        shell('rm -r swarm_vqtl_permute')
        
        # remove the temporary bed file (we already have gzipped copy)
        shell('rm data/moltraits-fastqtl.bed')
        
rule vqtlmapper__top_snp:
    """
    vqtl-mapper:    get top QTL per feature and correct for number of molecular
                    traits tested using Storey's FDR.
    """
    input:
        'vqtlmapper_permute.tsv.gz'
    output:
        'vqtlmapper_permute-top_qtl.tsv.gz'
    params:
        script_tophit=srcdir('scripts/general-min_line_per_feature.R'),
        script_fdr=srcdir('scripts/general-storey_fdr.R')
    run:
        # extract the most significant association for each mol. trait
        shell(
            'Rscript {params.script_tophit} '
                '--file {input} --id_column pheno_id --value_column P | '
            'Rscript {params.script_fdr} --column BETA --out_column q_beta | '
            'gzip -c > {output}'
        )
################################################################################

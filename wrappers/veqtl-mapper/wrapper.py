#!/usr/bin/env python3

from snakemake.shell import shell
import numpy as np


# the veqtl-mapper swarming mechanism is really silly...
# we need to count the number of lines in the intput file
with open(snakemake.input['pheno']) as f:
    n_features = sum(1 for line in f)
n_features -= 1 # drop header

# get the number of features we want to use per run
n_features_per_job = int(np.ceil(
    n_features / int(snakemake.wildcards['j_total'])))

# do we need the final job? -- could do this modulo
total_n = n_features_per_job * (int(snakemake.wildcards['j_total'])-1)
final_job = snakemake.wildcards['j_cur'] == snakemake.wildcards['j_total']
if total_n >= n_features and final_job:
    # we don't need the final job. example: 608 gene and 30 jobs
    shell('touch %s' % snakemake.output[0])
else:
    # start to build the command
    cmd = 'veqtl-mapper --vcf %s' % (snakemake.input['geno'])
    cmd = '%s --bed %s' % (cmd, snakemake.input['pheno'])
    cmd = '%s --genes %d' % (cmd, n_features_per_job) 
    cmd = '%s --job-number %s' % (cmd, snakemake.wildcards['j_cur']) 
    cmd = '%s --out %s' % (cmd, snakemake.output[0])

    params = dict(snakemake.params.items())
    
    window = params.get("window", "1000000") # cis window default = 1Mb
    cmd = '%s --window %s' % (cmd, window)
    
    if 'other_settings' in params:
        cmd = '%s %s' % (cmd, snakemake.params['other_settings'])
    
    if 'perm' in params:
        perm = '%s' % params['perm']
        if 'seed' in params:
            perm = '%s,%s' % (perm, params['seed'])
        cmd = '%s --perm %s' % (cmd, perm)
    
    shell(cmd)

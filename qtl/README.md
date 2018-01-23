Snakemake QTL pipeline
======================

Snakemake pipeline for molecular trait QTL mapping.

Overview
--------

The following tools are used for QTL mapping followed by the operations performed using each tool.

1. [QTLtools](https://qtltools.github.io/qtltools): map mean QTLs.
    
    A. Map nominal associations across all molecular traits (recording every variant tested). 
    
    B. Run 10,000 permutations for each molecular trait feature, controlling for FDR. 
    
    C. Map conditional QTLs using the thresholds identified for each molecular trait feature.

2. [veqtl-mapper](https://funpopgen.github.io/veqtl-mapper): map variance QTLs.
    
    A. Run 10,000 permutations for each molecular trait feature.
    
    * Papers: [Genetic interactions affecting human gene expression identified by variance association mapping](https://doi.org/10.7554/eLife.01381), [veqtl-mapper: variance association mapping for molecular phenotypes](https://doi.org/10.1093/bioinformatics/btx273)


Quickstart
----------

Quick demo shown below. 

Note that conda environment.yml contains a reproducible environment for most of the tools used in this repository except for QTLtools and veqtl-mapper, as they have not yet been built as a conda package. 

```bash
# make a directory for a demo
mkdir snkmk_qtl-demo
cd snkmk_qtl-demo

# set up the directory with the required data dir
mkdir data
cp ${SNK_REPO}/data/qtl/* data

# set up the config file inside the working dir
# change json file to fit the analysis
cp ${SNK_REPO}/qtl/config_analysis.json .

# now run snakemake in dryrun
snakemake --snakefile ${SNK_REPO}/qtl/Snakefile --configfile config_analysis.json --dryrun

# now run snakemake
snakemake --snakefile ${SNK_REPO}/qtl/Snakefile --configfile config_analysis.json --printshellcmds
```

The above jobs may take a little while, it would be much faster if we utilized a cluster to run the jobs. 

Below is a demo for LSF, you can easily tweak the config file for other systems.

```bash
# set up the config file inside the working dir
# change json files to fit the system
cp ${SNK_REPO}/config_cluster_lsf.json config_cluster.json

# now run snakemake
snakemake --snakefile ${SNK_REPO}/qtl/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster 'bsub -g {cluster.group} -M {cluster.memory} -R {cluster.resources} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}'

# alternatively use the cluster script (which will load settings from the config file)
snakemake --snakefile ${SNK_REPO}/qtl/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster ${SNK_REPO}/wrappers/cluster/lsf.py
```

Setup
-----

* The pipeline will run once the proper data frames are in the data directory inside the working directory. Use the demo data input files as a guide or see the QTLtools documentation. 
* Copy the config files in the working directory and change them as required. For instance, change the number of jobs to swarm various steps or change the number of permutations to 10000 (strongly suggested).


Working directory structure
---------------------------

* **Required input files**: in the data directory.
* **Final output data**: in the working directory.
* **Final summary plots**: in the plots directory.
* **Logs of all jobs**: in the logs directory.
* **Intermediate files**: currently all intermediate files are in the swarm* directories. These are deleted after the final output file is finished.


QTLtools
--------

* QTLtools can run with / without covariates and other PARAM. To change these across all scripts, edit the PARAM_QTLTOOLS parameter inside the config_analysis.json file. For instance, if you are using PEER residuals as the molecular trait input file, then then you would not want covariates (PARAM_QTLTOOLS=""). However, if you just used PEER to learn factors, you may want to include these factors as covariates in the QTLtools model. You can do that through the PARAM_QTLTOOLS parameter (PARAM_QTLTOOLS="--cov data/covariates.txt.gz"). Other PARAM could be to select only a subset of samples or molecular traits for analysis (e.g., PARAM_QTLTOOLS="--cov data/covariates.txt.gz --include-phenotypes data/moltraits-analyze.txt"). See the QTLtools [documentation](https://qtltools.github.io/qtltools) for various other PARAM.

* If you get an error, some helpful ideas:
    1. Are there cis variants around the molecular trait (e.g., you may have removed X chr variants from the vcf file and be trying to test for molecular traits on the X chr)? For nominal mapping, 0 cis variants will cause the molecular trait to be missing from the output file. For permutation mapping, the molecular trait will be included in the output but with NA values.
    
    2. Re-run a failed job in an interactive shell to see the exact failure. You can find the exact commands used by snakemake by looking at the logs in the /logs/cluster directory. 
    
    3. Is there something strange about the failed molecular trait readout (e.g, in the case of a gene, is it even expressed)?
    
    4. Finally, QTLtools bins jobs by breaking the chromsome into regions. If too many jobs are submitted and many phenotypes are excluded, you may find issues where a job iteration is linked to a region with no phenotypes. In this case, QTLtools will fail and the snakemake pipeline will not complete.


veqtl-mapper
------------

* Currently, veqtl-mapper does not support covariates. Therefore, it may be best to run veqtl-mapper on PEER residuals. 


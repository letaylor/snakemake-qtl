Snakemake PEER pipeline
=======================

Snakemake pipeline for latent factor analysis using PEER.


Overview
--------

1. [PEER](https://github.com/PMBio/peer).
    
    A. Run PEER using many different numbers of hidden factors. 
    
    B. [Optional] Optimize PEER for QTLs by mapping QTLs at each factor run. 

Note: When running this pipeline on a cluster, you may have to run snakemake multiple times or adjust the --latency-wait parameter. This is because some of the jobs take quite some time.


Quickstart
----------

Quick demo shown below. This assumes you ran the commands in the parent README file.

```bash
# make a directory for a demo
mkdir snkmk_peer-demo
cd snkmk_peer-demo

# set up the directory with the required data dir
mkdir data
cp ${SNK_REPO}/data/qtl/* data
gunzip -c data/covariates.txt.gz | sed s/' '/'\t'/g | gzip -c > data/covariates.tsv.gz 

# set up the config file inside the working dir
# change json file to fit the analysis
cp ${SNK_REPO}/peer/config_analysis.json .

# now run snakemake in dryrun
snakemake --snakefile ${SNK_REPO}/peer/Snakefile --configfile config_analysis.json --dryrun

# now run snakemake
snakemake --snakefile ${SNK_REPO}/peer/Snakefile --configfile config_analysis.json --printshellcmds
```

The above jobs may take a little while, it would be much faster if we utilized a cluster to run the jobs. 

Below is a demo for LSF, you can easily tweak the config file for other systems.

```bash
# set up the config file inside the working dir
# change json files to fit the system
cp ${SNK_REPO}/config_cluster_lsf.json config_cluster.json

# now run snakemake
snakemake --snakefile ${SNK_REPO}/peer/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster 'bsub -g {cluster.group} -M {cluster.memory} -R {cluster.resources} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}'

# alternatively use the cluster script (which will load settings from the config file)
snakemake --snakefile ${SNK_REPO}/peer/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster ${SNK_REPO}/wrappers/cluster/lsf.py

# finally, if we have optimized the number of PEER factors for QTL mapping, clear the logs for each PEER mapping
snakemake qtltools__in_subdir_clean_logs --snakefile ${SNK_REPO}/peer/Snakefile --configfile config_analysis.json --printshellcmd
```

Setup
-----

* The pipeline will run once the proper data frames are in the data directory inside the working directory. Use the demo data input files as a guide.
* Copy the config files in the working directory and change them as required. For instance, **change the number of iterations to 1000** (strongly suggested).


Working directory structure
---------------------------

* **Required input files**: in the data directory.
* **Final output data**: in the working directory and the peer_runs directory.
* **Final summary plots**: in the plots directory.
* **Logs of all jobs**: in the logs directory.


PEER
----

* PEER can run with / without covariates and other options. To change these across all scripts, edit the PARAM_PEER parameter inside the config_analysis.json file. If you do not wish to include covariates and would like to keep the defaults, then PARAM_PEER="". A list of options can be found by looking at the --help of scripts/peer-run.R.




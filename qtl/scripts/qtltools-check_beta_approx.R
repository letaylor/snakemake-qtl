#!/usr/bin/env Rscript

suppressMessages(library(data.table))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 2) 
    stop("Incorrect number of arguments, usage> Rscript qtltools-check_beta_approx.R INPUT_FILE OUTPUT"))

opt_input  = args[1];
opt_output = args[2];

d = data.frame(fread(paste("zcat", opt_input), sep="\t", header=T, 
    stringsAsFactors=F))


png(file=paste(opt_output, ".png", sep = ""), height=700, width=700, 
    units='px', res=120, type="cairo")
    
    plot(d$p_empirical, d$p_adjust_beta_dist, 
        xlab="Empirical p-value", 
        ylab="Beta approximation p-value", 
        main="")
    abline(0, 1, col="red")
    
dev.off()

png(file=paste(opt_output, "-neglog10.png", sep = ""), height=700, width=700, 
    units='px', res=120, type="cairo")
    
    plot(-log10(d$p_empirical), -log10(d$p_adjust_beta_dist), 
        xlab="Empirical p-value (-log10)", 
        ylab="Beta approximation p-value (-log10)", 
        main="")
    abline(0, 1, col="red")
    
dev.off()
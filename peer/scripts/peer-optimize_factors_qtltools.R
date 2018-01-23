#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

######################## Required Packages #####################################
suppressMessages(library(ggplot2))
################################################################################

######################## Functions #############################################
run = function() {
  suppressPackageStartupMessages(library(optparse))
  
  optionList = list(
    make_option(c("-f", "--file"), 
      type="character", 
      help=paste0(
        "Input tab delimited file with the following columns: ",
        "(peer_factors, fdr, n) or (peer_factors, p_cut, n)."
      )
    ),
    
    make_option(c("-y", "--y_axis_label"), 
      type="character", default="Number of molecular traits with a QTL", 
      help="Label for y axis. [default %default]"
    ),
    
    make_option(c("-o", "--out_file"), 
      type="character", default="peer-optimize_factors_QTLtools", 
      help="Name of output file (will be in pdf format)"
    )
  )
  
  parser = OptionParser(
      usage='%prog', 
      option_list=optionList,
      description = paste0(
          'Plots summary of QTLs per PEER factor'
      )
  )
  
  # a hack to fix a bug in optparse that won't let you use positional args
  # if you also have non-boolean optional args:
  getOptionStrings = function(parserObj) {
      optionStrings = character()
      for(item in parserObj@options) {
          optionStrings = append(optionStrings, 
              c(item@short_flag, item@long_flag))
      }
      optionStrings
  }
  optStrings = getOptionStrings(parser)
  arguments = parse_args(parser, positional_arguments=TRUE)
  
  
  # read in the matrix
  f=arguments$options$file
  if (f == 'stdin') {
    dat = read.table(file=file('stdin'), sep='\t', header=TRUE, 
      stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    dat = read.table(file=f, header=TRUE, sep='\t')
  }
  
  color_lab="FDR"
  if ("p_cut" %in% colnames(dat)) {
    dat$fdr=dat$p_cut
    color_lab="p-value cutoff"
  }
  
  base=arguments$options$out_file
  y_lab=arguments$options$y_axis_label
  pdf(file=paste(base, ".pdf", sep = ""), height=4.5, width=5)
    plt = ggplot(dat, aes(x=peer_factors, y=n, color=as.factor(fdr))) 
    plt = plt + theme_bw() + scale_color_brewer(palette="Dark2")
    plt = plt + geom_point(alpha=0.4) + geom_line(alpha=0.7)
    plt = plt + labs(x="PEER factors", y=y_lab, color=color_lab)
    print(plt)
  dev.off()
  
  return(0)
}
################################################################################


run()

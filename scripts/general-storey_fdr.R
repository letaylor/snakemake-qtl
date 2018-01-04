#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList = list(
    make_option(c('-f', '--file'), 
        type='character',
        default='stdin',
        help=paste0(
            'Input file name or stdin if reading from a pipe. ',
            'Delimiter is "\\t". [default %default]'
        )
    ), 
        
    make_option(c('-c', '--column'), 
        type='character', 
        help='Column to use for input values.'
    ), 
        
    make_option(c('-oc', '--out_column'), 
        type='character',
        default='',
        help='Output column name. [default q_<column>]'
    )
)

parser = OptionParser(
    usage='%prog -f file -c column -oc output_column_name', 
    option_list=optionList,
    description=paste0(
        'Performs Storey\'s FDR procedure on a dataframe. ',
        'Output returned to stdout.'
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


arguments = parse_args(parser, positional_arguments=TRUE)
################################################################################


######################## Required Packages #####################################
suppressPackageStartupMessages(library(qvalue))
################################################################################


######################## Read Data & Manipulate ################################
# read in the input data
f=arguments$options$file
if (f == 'stdin') {
    dat = read.table(file=file('stdin'), sep="\t", header=TRUE, 
        stringsAsFactors=FALSE)
} else {
    #dat = read.table(file=f, header=T, sep='\t', stringsAsFactors=FALSE)
    suppressMessages(library(data.table))
    dat = data.frame(fread(paste("zcat", f), sep="\t", header=TRUE, 
        stringsAsFactors=FALSE, showProgress=FALSE))
}

# set up in and out columns
in_col = arguments$options$column
if (arguments$options$out_column == '') {
    out_col = paste('q', in_col, sep='__')
} else {
    out_col = arguments$options$out_column
}

dat[[out_col]] = qvalue(dat[[in_col]])$qvalues

# output to a file
write.table(dat,
    file=stdout(),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t")
################################################################################


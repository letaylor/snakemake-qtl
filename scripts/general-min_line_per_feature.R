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
        
    make_option(c('-ic', '--id_column'), 
        type='character', 
        help='Column to use for ids.'
    ), 
        
    make_option(c('-vc', '--value_column'), 
        type='character',
        help='Column to use for values'
    )
)

parser = OptionParser(
    usage='%prog -f file -ic id_column -vc value_column', 
    option_list=optionList,
    description = paste0(
        'Finds the minimum line for each id based on value. ',
        'For ties returns the first instance.'
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
suppressPackageStartupMessages(library(plyr))
################################################################################


######################## Read Data & Manipulate ################################
# read in the input data
f=arguments$options$file
if (f == 'stdin') {
    dat = read.table(file=file('stdin'), sep='\t', header=TRUE, 
        stringsAsFactors=FALSE)
} else {
    suppressMessages(library(data.table))
    #dat = read.table(file=f, header=T, sep='\t')
    dat = data.frame(fread(paste("zcat", f), sep="\t", header=TRUE, 
        stringsAsFactors=FALSE, showProgress=FALSE))
}

# set up in and out columns
id_col = arguments$options$id_column
val_col = arguments$options$value_column

# make temp variables
dat[["val_col"]] = dat[[val_col]] 

# get the min value across the id col
dat_out = plyr::ddply(dat, c(id_col), function (df) {
    d=subset(df, val_col == min(df$val_col, na.rm=TRUE));
    return(d[1,]) # return first row (in case of ties)
})

# delete temp variable
dat_out[["val_col"]] = NULL 

# output to a file
write.table(dat_out,
    file=stdout(),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t")
################################################################################


print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(stringr)
    library(argparse)
})


#COPIED from Signac::BinarizeCounts
makeBinary <- function(object) {
    if (inherits(x = object, what = "CsparseMatrix")) {
        slot(object = object, name = "x") <- rep.int(x = 1, times = length(x = slot(object = object, name = "x")))
    } else {
        object[object > 1] <- 1
    }
    return(object)
}


parser <- ArgumentParser(description='Subset and reorder matrices')
parser$add_argument('mat_file', metavar='mat_file', help='feature x cells file')
parser$add_argument('meta_file', metavar='meta_file', help='cells x metadata file')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
meta_file <- args$meta_file
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time('Load files')
mat <- readRDS(mat_file)
dim(mat)
mat[1:3,1:3]

meta <- readRDS(meta_file)
dim(meta)
meta[1:3,1:3]


print_time('Checking cells')
if(!any(rownames(meta) %in% rownames(mat))){
    cat('Meta cell names did not match matrix cell names; trying ArchR style.\n')
    
    if(all(grepl('-1$',rownames(meta)))){
        cellNames <- sub('#','_',rownames(mat))
    } else {
        split <- str_split_fixed(rownames(mat),'#',2)
        cellNames <- paste(sep='_', split[,1],str_split_fixed(split[,2],'-',2)[,1])
    }
    rownames(mat) <- cellNames
}
if(!identical(sort(rownames(meta)),sort(rownames(mat)))) stop('matrix and meta rownames need to be the same.')


print_time('Reorder cells')
mat <- mat[rownames(meta),]
dim(mat)
mat[1:3,1:3]
if(!identical(rownames(meta),rownames(mat))) stop('Cells not identical')


print_time('Saving')
saveRDS(mat,paste(sep="",outPrefix,'.rds'))


print_time('Done.')


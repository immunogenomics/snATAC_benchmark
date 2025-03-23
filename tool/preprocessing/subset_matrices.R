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
parser$add_argument('--rmZeroFeat', action='store_true', help='if given, remove features with zero counts')
parser$add_argument('--binarize', action='store_true', help='if given, output binary matrix')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
meta_file <- args$meta_file
rmZeroFeat <- args$rmZeroFeat
binarize <- args$binarize
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
if(!any(rownames(meta) %in% colnames(mat))){
    cat('Meta cell names did not match matrix cell names; trying ArchR style.\n')
    split1 <- str_split_fixed(colnames(mat),'#',2)
    split2 <- str_split_fixed(split1[,2],'-',2)
    colnames(mat) <- paste(sep='',split1[,1],'_',split2[,1])
}
if(!all(rownames(meta) %in% colnames(mat))) stop('All meta rownames need to be in matrix rownames.')


print_time('Subset cells')
mat <- mat[,rownames(meta)]
dim(mat)
if(!identical(rownames(meta),colnames(mat))) stop('Cells not identical')


if(rmZeroFeat){
    print_time('Remove all zero features')
    
    rs <- rowSums(mat)
    print(head(table(rs)))
    mat <- mat[names(rs[rs!=0]),]
    print(dim(mat))
    rs <- rowSums(mat)
    print(head(table(rs)))
    
    outPrefix <- paste(sep="",outPrefix,'_nonzero')
}


if(binarize){
    print_time('Binarize')

    cat(paste('Original max:',max(mat),'\n'))
    mat <- makeBinary(mat)
    if(max(mat)!=1) stop('Binarization failed')
    cat(paste('New max:',max(mat),'\n'))
    
    outPrefix <- paste(sep="",outPrefix,'_binary')
}


print_time('Saving')
saveRDS(mat,paste(sep="",outPrefix,'.rds'))


print_time('Done.')


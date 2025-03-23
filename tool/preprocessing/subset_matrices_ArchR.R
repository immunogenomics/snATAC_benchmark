print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
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
parser$add_argument('ArchR_sumExp_file', metavar='ArchR_sumExp_file', help='feature x cells summarized experiment file from ArchR; cells: sample#CB-1')
parser$add_argument('assay', metavar='assay', help='assay to extract - e.g., GeneScoreMatrix')
parser$add_argument('meta_file', metavar='meta_file', help='cells x metadata file')
parser$add_argument('--delim', metavar='delim', default='_', help='cell name delimiter between sample and cell barcode; default _')
parser$add_argument('--rmZeroFeat', action='store_true', help='if given, remove features with zero counts')
parser$add_argument('--binarize', action='store_true', help='if given, output binary matrix')
parser$add_argument('--NOrm1', action='store_true', help='if given, do NOT remove the "-1" at the end of the cell barcode')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

ArchR_sumExp_file <- args$ArchR_sumExp_file
assay <- args$assay
meta_file <- args$meta_file
delim <- args$delim
rmZeroFeat <- args$rmZeroFeat
binarize <- args$binarize
NOrm1 <- args$NOrm1
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(ArchR_sumExp_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time('Load sumExp')
ArchR_sumExp <- readRDS(ArchR_sumExp_file)
dim(ArchR_sumExp)
ArchR_sumExp


print_time('Process sumExp')
if(!(assay %in% names(assays(ArchR_sumExp)))) stop('This assay does not exist in this Summarized Experiment.')
mat <- assays(ArchR_sumExp)[[assay]]

df <- as.data.frame(rowData(ArchR_sumExp),stringsAsFactors=FALSE)
if('name' %in% colnames(df)){
    rownames(mat) <- df$name
} else if('end' %in% colnames(df)){
    df$feat <- paste(sep='',df$seqnames,':',df$start,'-',df$end)
    rownames(mat) <- df$feat
} else if('start' %in% colnames(df)){
    tileSize <- df[2,'start']-df[1,'start']
    df$feat <- paste(sep='',df$seqnames,':',df$start,+tileSize)
    rownames(mat) <- df$feat
}

mat <- mat[,sort(colnames(mat))]
if(NOrm1){
    cellNames <- sub('#',delim,colnames(mat))
} else {
    split <- str_split_fixed(colnames(mat),'#',2)
    cellNames <- paste(sep=delim, split[,1],str_split_fixed(split[,2],'-',2)[,1])
}
colnames(mat) <- cellNames

dim(mat)
mat[1:3,1:3]


print_time('Checking cells')
meta <- readRDS(meta_file)
dim(meta)
meta[1:3,1:3]

if(!all(rownames(meta) %in% colnames(mat))) stop('All meta rownames need to be in matrix rownames.')


print_time('Subset cells')
mat <- mat[,rownames(meta)]
dim(mat)
mat[1:3,1:3]
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


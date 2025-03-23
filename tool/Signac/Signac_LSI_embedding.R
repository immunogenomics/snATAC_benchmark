print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Seurat)
    library(Signac)
    library(GenomicRanges)
    library(argparse)
})

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}

#COPIED from Signac::BinarizeCounts
makeBinary <- function(object) {
    if (inherits(x = object, what = "CsparseMatrix")) {
        slot(object = object, name = "x") <- rep.int(x = 1, times = length(x = slot(object = object, name = "x")))
    } else {
        object[object > 1] <- 1
    }
    return(object)
}


parser <- ArgumentParser(description='Signac LSI Embedding')
parser$add_argument('mat_file', metavar='mat_file', help='feat x cells file')
parser$add_argument('--binarize', action='store_true', help='if given, binarize matrix')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('--sample_col', metavar='sample_col', default='sample', help='sample column; default sample')
parser$add_argument('--rmSample_lessN', metavar='rmSample_lessN', type='integer', default=-1, help='remove samples with less than X cells; use -1 to do nothing; default -1')
parser$add_argument('--max_dim', metavar='max_dim', type='integer', default=30, help='max dimensions; default 30')
parser$add_argument('--NOrmDim1', action='store_true', help='if given, do not remove first dimension - often correlated with fragment count')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
binarize <- args$binarize
seed <- args$seed
meta_file <- args$meta_file
sample_col <- args$sample_col
rmSample_lessN <- args$rmSample_lessN
max_dim <- args$max_dim
NOrmDim1 <- args$NOrmDim1
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_LSI_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)

print_time('Load mat file')
mat <- readRDS(mat_file)
max(mat)
dim(mat)
mat[1:3,1:3]
if(binarize){
    mat <- makeBinary(mat)
    cat(paste('Binarized. New max:',max(mat),'\n'))
}

print_time('Load meta file')
meta <- readRDS(meta_file)
if(!identical(sort(rownames(meta)),sort(colnames(mat)))){
    stop('cells must be the same between meta and mat')
} else if(!identical(rownames(meta),colnames(mat))){
    meta <- meta[colnames(mat),]
    if(!identical(rownames(meta),colnames(mat))) stop('cells are not the same.')
} else {
    #do nothing since cellnames are identical
    doNothing = TRUE
}
dim(meta)
head(meta,n=3)
if(!(sample_col %in% colnames(meta))) stop('sample_col must be in meta file')

if(rmSample_lessN!=-1){
    cat(paste('Removing samples with less than',rmSample_lessN,'cells.\n'))
    sample_count <- as.data.frame(table(meta[,sample_col]),stringsAsFactors=FALSE)
    sample_toKeep <- sample_count[which(sample_count$Freq>=rmSample_lessN),'Var1']
    meta <- meta[which(meta[,sample_col] %in% sample_toKeep),]
    mat <- mat[,rownames(meta)]
    if(!identical(rownames(meta),colnames(mat))) stop('cell names not identical after subsetting.')
    cat(paste('Using',nrow(meta),'cells.\n'))
}


print_time('Create Seurat Object')
obj <- CreateSeuratObject(counts = mat,assay = "peaks",meta.data = meta)
obj

print_time('LSI processing')
obj <- FindTopFeatures(obj)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj,n = max_dim)
obj


print_time('Saving')
emb <- Embeddings(obj, reduction = "lsi")
if(!NOrmDim1){
    emb <- emb[,2:max_dim]
}
dim(emb)
emb[1:3,1:3]
saveRDS(emb,paste(sep='',outPrefix,'_LSI_embeddings.rds'))

saveRDS(obj,paste(sep='',outPrefix,'_LSI_object.rds'))


print_time("Done")

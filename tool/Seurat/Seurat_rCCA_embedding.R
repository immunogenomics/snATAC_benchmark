print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Seurat)
    library(scales)
    library(ggplot2)
    library(patchwork)
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

parser <- ArgumentParser(description='Seurat rCCA Embedding')
parser$add_argument('mat_file', metavar='mat_file', help='feat x cells file')
parser$add_argument('--is_norm', action='store_true', help='if given, mat file is already normalized')
parser$add_argument('--do_not_log', action='store_true', help='if given, do not log normalized mat file with log1p - only used if is_norm is TRUE')
parser$add_argument('--varFeat_method', metavar='varFeat_method', default='vst', choices=c('vst','mvp','disp'), help='FindVariableFeatures selection.method - options: vst, mvp, disp; default vst')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('--sample_col', metavar='sample_col', default='sample', help='sample column; default sample')
parser$add_argument('--rmSample_lessN', metavar='rmSample_lessN', type='integer', default=-1, help='remove samples with less than X cells; use -1 to do nothing; default -1')
parser$add_argument('--max_dim', metavar='max_dim', default=30, type='integer', help='max dimensions; default: 30')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
is_norm <- args$is_norm
do_not_log <- args$do_not_log
varFeat_method <- args$varFeat_method
meta_file <- args$meta_file
sample_col <- args$sample_col
rmSample_lessN <- args$rmSample_lessN
max_dim <- args$max_dim
seed <- args$seed
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_rCCA_args.txt")
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

print_time('Load meta file')
meta <- readRDS(meta_file)
if(!identical(sort(rownames(meta)),sort(colnames(mat)))){
    stop('cells must be the same between meta and mat')
} else if(!identical(rownames(meta),colnames(mat))){
    meta <- meta[colnames(mat),]
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
if(!is_norm){
    obj <- CreateSeuratObject(counts = mat, meta.data=meta)
} else {
    if(do_not_log){
        assay5 <- CreateAssay5Object(data=mat)
    } else {
        assay5 <- CreateAssay5Object(data=log1p(mat))
    }
    obj <- CreateSeuratObject(assay5,meta.data=meta)
}
obj


print_time('Split by sample')
ll <- obj[[sample_col]]
samp_vec <- ll[,1]
names(samp_vec) <- rownames(ll)

obj[["RNA"]] <- split(obj[["RNA"]], f = samp_vec)


print_time("Unintegrated processing")
if(!is_norm){
    obj <- NormalizeData(obj)
}
obj <- FindVariableFeatures(obj,selection.method=varFeat_method)
obj <- ScaleData(obj)
obj <- RunPCA(obj, seed.use=seed, npcs = max_dim)


print_time("Integrated processing")
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE, dims=1:max_dim
)


print_time("Saving")
emb <- Embeddings(obj, reduction = "integrated.cca")
dim(emb)
emb[1:3,1:3]
saveRDS(emb,paste(sep='',outPrefix,'_rCCA_embeddings.rds'))

saveRDS(obj,paste(sep='',outPrefix,'_rCCA_object.rds'))


print_time("Done")

print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(harmony)
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


parser <- ArgumentParser(description='Harmony correction of existing embedding')
parser$add_argument('embMat_file', metavar='embMat_file', help='feature x dim file')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('harmony_covars', metavar='harmony_covars', help='meta columns to use in harmony')
parser$add_argument('--max_dim', metavar='--max_dim', type='integer', help='max dimensions; default use all provided')
parser$add_argument('--harmony_sigma', metavar='harmony_sigma', type='double', default = 0.1, help='width of soft kmeans clusters in Harmony; bigger = more diffuse; default 0.1')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

embMat_file <- args$embMat_file
seed <- args$seed
meta_file <- args$meta_file
harmony_covars <- args$harmony_covars
max_dim <- args$max_dim
harmony_sigma <- args$harmony_sigma
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(embMat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_Harmony_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)


print_time('Load embeddings')
embMat <- readRDS(embMat_file)
if(!is.null(max_dim)){
    if(max_dim>ncol(embMat)){
        stop(paste('Asked for',max_dim,'dimensions; only have',ncol(embMat)))
    } else if(max_dim<ncol(embMat)){
        cat(paste('Subsetting from',ncol(embMat),'dimensions to the asked for',max_dim,'\n'))
        embMat <- embMat[,1:max_dim]
    } #if equal, do nothing!
}
dim(embMat)
embMat[1:3,1:3]


print_time('Load meta file')
meta <- readRDS(meta_file)
harmony_col_toUse <- parse_columns(harmony_covars)
if(!all(harmony_col_toUse %in% colnames(meta))) stop('not all harmony covariates in meta file')


print_time('Harmony')
set.seed(seed)
harmonyObj <- RunHarmony(embMat, meta, harmony_col_toUse, sigma = harmony_sigma,
                         return_object=TRUE)

harmony_mat <- as.matrix(t(harmonyObj$Z_corr))
rownames(harmony_mat) <- rownames(embMat)
colnames(harmony_mat) <- paste(sep='','h',colnames(embMat))
saveRDS(harmony_mat,paste(sep='',outPrefix,'_Harmony_embeddings.rds'))


print_time('Done.')


print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(irlba)
    library(argparse)
})


parser <- ArgumentParser(description='PCA of a feature selected file')
parser$add_argument('featSel_file', metavar='featSel_file', help='Feature selection file (cells x feature) - assumes it is center/scale already')
parser$add_argument('--centerScale', action='store_true', help='If given, center/scale')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('--max_dim', metavar='max_dim', type='integer', default = 30, help='number of PCs to calculate; default 30')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

featSel_file <- args$featSel_file
centerScale <- args$centerScale
seed <- args$seed
max_dim <- args$max_dim
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(featSel_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_PCA_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)


print_time('Load feature selection file')
featSel_mat <- readRDS(featSel_file)
if(centerScale){
    featSel_mat <- scale(featSel_mat,center=TRUE)
}
featSel_mat[1:3,1:3]
cat(paste("first column: mean:",mean(featSel_mat[,1]),"var: ",var(featSel_mat[,1]),"\n"))
summary(featSel_mat[,1])


print_time("PCA")
set.seed(seed)
pcaRes <- prcomp_irlba(featSel_mat, n = max_dim, center = FALSE, scale = FALSE) #matrix already has center/scale
saveRDS(pcaRes,paste(sep='',outPrefix,'_PCA_res.rds'))

pcaMat <- pcaRes$x
rownames(pcaMat) <- rownames(featSel_mat)
dim(pcaMat)
pcaMat[1:3,1:3]
saveRDS(pcaMat,paste(sep='',outPrefix,'_PCA_embeddings.rds'))


print_time('Done.')


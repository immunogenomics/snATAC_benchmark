print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(RANN)
    library(argparse)
})


makeNNgraph <- function (data.use, k.param = 30, nn.eps = 0)
{
    my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    
    nn_mat <- t(apply(nn.ranked,1,function(x){sort(rownames(data.use[x,]))}))
    rownames(nn_mat) <- rownames(data.use)
    colnames(nn_mat) <- paste(sep='','NN',1:k.param)
    
    return(nn_mat)
}


parser <- ArgumentParser(description='NxK KNN matrix')
parser$add_argument('emb_file', metavar='emb_file', help='cell x dim file')
parser$add_argument('--knn', metavar='knn', type='integer', action='append', help='number of nn for NN graph - can input multiple knn values')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

emb_file <- args$emb_file
knn <- args$knn
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(emb_file))) stop("Input file doesn't exist.")
if(length(knn)==0) stop('Must have at least one knn value.')

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_NNgraph_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time('Load embedding file')
embMat <- readRDS(emb_file)
dim(embMat)
embMat[1:3,1:3]


print_time('KNN graph')
for(kval in knn){
    NNgraph <- makeNNgraph(embMat,k.param=kval)
    saveRDS(NNgraph, paste(sep='',outPrefix,'_NN',kval,'graph.rds'))
}

print_time('Done.')


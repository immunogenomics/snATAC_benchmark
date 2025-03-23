print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Seurat)
    library(parallel)
    library(argparse)
})

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0)
{
    my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard",
        eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn_res) <- row.names(data.use)
    colnames(snn_res) <- row.names(data.use)
    return(snn_res)
}
environment(BuildSNNSeurat) <- asNamespace("Seurat")

getClusters <- function(adt_harmony_res, resolution_list){
    snn_ref <- BuildSNNSeurat(adt_harmony_res, nn.eps = 0) 

    ids_ref <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
        Seurat:::RunModularityClustering(SNN = snn_ref, modularity = 1,
            resolution = res_use, algorithm = 1, n.start = 20,
            n.iter = 20, random.seed = 100, print.output = FALSE,
            temp.file.location = NULL, edge.file.name = NULL)
    }, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list)))) 

    ids_ref <- as.data.frame(ids_ref)
    rm(snn_ref)
    colnames(ids_ref) <- sprintf("res_%.2f", resolution_list)
    rownames(ids_ref) <- rownames(adt_harmony_res)

    return(ids_ref)
}

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}


parser <- ArgumentParser(description='Make clusters from embedding matrix')
parser$add_argument('embMat_file', metavar='embMat_file', help='feature x dim file')
parser$add_argument('--max_dim', metavar='--max_dim', type='integer', help='max dimensions; default use all provided')
parser$add_argument('--cluster_res', metavar='cluster_res', default='0.2,0.4,0.6,0.8,1.0', help='cluster resolutions, separated by comma; default 0.2,0.4,0.6,0.8,1.0')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

embMat_file <- args$embMat_file
max_dim <- args$max_dim
cluster_res <- args$cluster_res
seed <- args$seed
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(embMat_file))) stop("Input file doesn't exist.")

res_list <- as.numeric(parse_columns(cluster_res))
if(any(is.na(res_list))) stop('non-numeric cluster resolutions given!')
paste('Cluster Resolutions:',paste(collapse=', ',res_list))

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_clusters_args.txt")
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


print_time('Clustering')
set.seed(seed)
cluster_res <- getClusters(embMat,res_list)
cluster_df <- data.frame(apply(cluster_res, 2, as.character),stringsAsFactors=FALSE)
rownames(cluster_df) <- rownames(embMat)
saveRDS(cluster_df,paste(sep='',outPrefix,'_dim',ncol(embMat),'_cluster.rds'))


print_time('Done.')


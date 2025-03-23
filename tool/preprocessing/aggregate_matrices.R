print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    #library(Matrix.utils)
    library(stringr)
    library(argparse)
})

aggMat <- function(thisMat,thisMeta,thisCTcol,fun){

    theseCT <- sort(unique(thisMeta[,thisCTcol]))
    for(idx in 1:length(theseCT)){
        thisCT <- theseCT[idx]
        these_cells <- rownames(thisMeta[which(thisMeta[,thisCTcol]==thisCT),])
        this_subset <- thisMat[,these_cells]

        if(length(these_cells)==1){
            thisAgg <- this_subset
        } else {
            thisAgg <- as.matrix(apply(this_subset,1,fun))
        }

        if(idx==1){
            catAgg <- thisAgg
        } else {
            catAgg <- cbind(catAgg,thisAgg)
        }
        colnames(catAgg)[ncol(catAgg)] <- thisCT
    }
    
    return(as(catAgg,'dgCMatrix'))

}

agg_merge <- function(half1, half2){
    return(rbind(half1,half2))
}       
    
agg_byHalf <- function(thisMat,thisMeta,thisCTcol,fun,num_rows_max=24000){
    if(nrow(thisMat)<=num_rows_max){
        return(aggMat(thisMat,thisMeta,thisCTcol,fun))
    } else {
        mid_row_idx <- floor(nrow(thisMat)/2)
        a <- agg_byHalf(thisMat[1:mid_row_idx,],thisMeta,thisCTcol,fun,num_rows_max=num_rows_max)
        b <- agg_byHalf(thisMat[(mid_row_idx+1):nrow(thisMat),],thisMeta,thisCTcol,fun,num_rows_max=num_rows_max)
        agg_merge(a,b)
    }
}

checkMat <- function(thisMat){
    rs <- rowSums(thisMat)
    cat(paste(sep='','Number of all zero rows: ',length(rs[rs==0]),'\n'))

    cs <- colSums(thisMat)
    cat(paste(sep='','Number of all zero cell (types): ',length(cs[cs==0]),'\n'))
}


parser <- ArgumentParser(description='Mean-aggregate matrix by meta column')
parser$add_argument('mat_file', metavar='mat_file', help='matrix file')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('col_name', metavar='col_name', help='cell type column')
parser$add_argument('--nrow_forHalf', metavar='nrow_forHalf', type='integer', default=200000, help='max number of rows to give to aggregate matrix; default 200000')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
meta_file <- args$meta_file
col_name <- args$col_name
nrow_forHalf <- args$nrow_forHalf
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_meanAggr_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load matrix")
mat <- readRDS(mat_file)
dim(mat)
max(mat)
mat[1:3,1:3]
checkMat(mat)

print_time("Load meta")
meta <- readRDS(meta_file)
if(!(col_name %in% colnames(meta))) stop('given column name not in meta_file')
if(!all(colnames(mat) %in% rownames(meta))) stop('not all matrix cells are in meta')
meta <- meta[colnames(mat),]
if(!identical(colnames(mat),rownames(meta))) stop('cells not identical between matrix and meta')


print_time('Mean-aggregate')
mat_aggr <- agg_byHalf(mat,meta,col_name,mean,num_rows_max=nrow_forHalf)
saveRDS(mat_aggr,paste(sep='',outPrefix,'_meanAggr.rds'))
checkMat(mat_aggr)

print_time('Done.')


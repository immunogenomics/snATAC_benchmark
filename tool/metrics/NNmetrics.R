print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(stringr)
    library(stringi)
    library(philentropy)
    library(argparse)
})


split_sample <- function(cb,delimiter='_'){
    return(stri_reverse(str_split_fixed(stri_reverse(cb),'_',2)[,2]))
}


NN_metrics <- function(gene_KNN,peak_KNN,allProb,uniqSamples,delimiter='_',removeItself=TRUE,mType='NxN',k.param = 100){
    if(!identical(rownames(gene_KNN),rownames(peak_KNN))) stop('rownames must be identical')
    if(!(mType %in% c('NxK','NxN'))) stop('matrix type needs to be NxK or NxN')
    
    
    metrics_df <- data.frame('sample'=character(),'mKNN'=integer(),'lsKLD'=numeric(),stringsAsFactors=FALSE)
    

    for(idx in 1:nrow(gene_KNN)){
        thisCell <- rownames(gene_KNN)[idx]
        thisSample <- split_sample(thisCell,delimiter=delimiter)
        
        if(mType=='NxN'){
            gnn <- gene_KNN[idx,]
            gnn <- names(gnn[gnn==1])

        } else {
            gnn <- unname(gene_KNN[idx,])
        }
        if(length(gnn)!=k.param){
            cat(paste('WARNING: gene',idx,'does not have',k.param,'NN'))
        }
        if(removeItself){
            gnn <- gnn[gnn!=thisCell]
        }
        samples_thisCellKNN_gene <- split_sample(gnn,delimiter=delimiter)
        bySample_count_thisCellKNN_gene <- table(factor(samples_thisCellKNN_gene,levels=uniqSamples))
        bySample_count_thisCellKNN_gene_psAllProb <- bySample_count_thisCellKNN_gene+allProb
        bySample_prob_thisCellKNN_gene <- bySample_count_thisCellKNN_gene_psAllProb/sum(bySample_count_thisCellKNN_gene_psAllProb)
        if(!all.equal(sum(bySample_prob_thisCellKNN_gene),1)) stop('gene probabilities must sum to 1')
        if(any(bySample_prob_thisCellKNN_gene==0)) stop('some gene probabilities are 0')
        
        if(mType=='NxN'){
            pnn <- peak_KNN[idx,]
            pnn <- names(pnn[pnn==1])

        } else {
            pnn <- unname(peak_KNN[idx,])
        }
        if(length(pnn)!=k.param){
            cat(paste('WARNING: peak',idx,'does not have',k.param,'NN'))
        }
        if(removeItself){
            pnn <- pnn[pnn!=thisCell]
        }
        samples_thisCellKNN_peak <- split_sample(pnn,delimiter=delimiter)
        bySample_count_thisCellKNN_peak <- table(factor(samples_thisCellKNN_peak,levels=uniqSamples))
        bySample_count_thisCellKNN_peak_psAllProb <- bySample_count_thisCellKNN_peak+allProb
        bySample_prob_thisCellKNN_peak <- bySample_count_thisCellKNN_peak_psAllProb/sum(bySample_count_thisCellKNN_peak_psAllProb)
        if(!all.equal(sum(bySample_prob_thisCellKNN_peak),1)) stop('peak probabilities must sum to 1')
        if(any(bySample_prob_thisCellKNN_peak==0)) stop('some peak probabilities are 0')
        

        mKNN <- length(gnn[which(gnn %in% pnn)])
        
        lsKLD <- suppressMessages(unname(philentropy::KL(rbind(bySample_prob_thisCellKNN_gene,
                                                               bySample_prob_thisCellKNN_peak))))
        

        metrics_df <- rbind(metrics_df,data.frame('sample'=thisSample,'mKNN'=mKNN,'lsKLD'=lsKLD,
                                             stringsAsFactors=FALSE,row.names=c(thisCell)))
    }
    
    return(metrics_df)

}



parser <- ArgumentParser(description='snATAC-seq benchmark NN-based metrics')
parser$add_argument('gene_knn_file', metavar='gene_knn_file', help='RNA KNN file; either NxN or NxK')
parser$add_argument('peak_knn_file', metavar='peak_knn_file', help='ATAC KNN file; either NxN or NxK')
parser$add_argument('--knn', metavar='knn', default=100, type='integer', help='number of NN; default 100')
parser$add_argument('--delim', metavar='delim', default='_', help='delimiter between sample and barcode in cell rownames - allows for multiple delimiters in sample; default _')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

gene_knn_file <- args$gene_knn_file
peak_knn_file <- args$peak_knn_file
knn <- args$knn
delim <- args$delim
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(gene_knn_file,peak_knn_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_NNmetrics_args.txt")

for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time('Load KNN files')
gKNN <- readRDS(gene_knn_file)
pKNN <- readRDS(peak_knn_file)

if(!identical(rownames(gKNN),rownames(pKNN))){
    gKNN <- gKNN[sort(rownames(gKNN)),]
    pKNN <- pKNN[sort(rownames(pKNN)),]

    if(!identical(rownames(gKNN),rownames(pKNN))){
        stop('ERROR: Not the same cells in the gene and peak KNN!')
    }
}

matType <- 'neither'
if(nrow(gKNN)==ncol(gKNN) & nrow(pKNN)==ncol(pKNN)){
    matType <- 'NxN'
} else if(ncol(gKNN)==knn & ncol(pKNN)==knn){
    matType <- 'NxK'
} else {
    stop('input KNN matrices must be NxN or NxK')
}


print_time('All Cells sample probabilities')
bySample_count_cells_all <- table(split_sample(rownames(gKNN),delimiter=delim))
bySample_prob_cells_all <- bySample_count_cells_all/sum(bySample_count_cells_all)
if(!all.equal(sum(bySample_prob_cells_all),1)) stop('all cells probabilities must sum to 1')


print_time('Per-cell sample probabilities')
ret <- NN_metrics(gKNN,pKNN,bySample_prob_cells_all,names(bySample_count_cells_all),
                  delimiter=delim,removeItself=TRUE,mType=matType,k.param=knn)
saveRDS(ret,paste(sep='',outPrefix,'_NNmetrics.rds'))

summary(ret$mKNN)
summary(ret$lsKLD)


print_time('Done.')


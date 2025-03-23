print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Seurat)
    library(symphony)
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

#COPIED from Seurat v3 TF.IDF
TF.IDF <- function (data, verbose = TRUE) {
    if (is.data.frame(x = data)) {
        data <- as.matrix(x = data)
    }
    if (!inherits(x = data, what = "dgCMatrix")) {
        data <- as(object = data, Class = "dgCMatrix")
    }
    if (verbose) {
        message("Performing TF-IDF normalization")
    }
    npeaks <- colSums(x = data)
    tf <- t(x = t(x = data)/npeaks)
    idf <- ncol(x = data)/rowSums(x = data)
    norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
    norm.data[which(x = is.na(x = norm.data))] <- 0
    rownames(norm.data) <- rownames(data) #added
    return(norm.data)
}

scale_merge <- function(half1, half2){
    return(cbind(half1,half2))
}

scale_byHalf <- function(mat,num_cols_max=24000,cntr=TRUE,scl=TRUE){
    if(ncol(mat)<=num_cols_max){
        return(scale(mat,center=cntr,scale=scl))
    } else {
        mid_col_idx <- floor(ncol(mat)/2)
        a <- scale_byHalf(mat[,1:mid_col_idx],num_cols_max=num_cols_max,cntr=cntr,scl=scl)
        b <- scale_byHalf(mat[,(mid_col_idx+1):ncol(mat)],num_cols_max=num_cols_max,cntr=cntr,scl=scl)
        scale_merge(a,b)
    }
}


parser <- ArgumentParser(description='minially accessible, norm, and variable peaks using symphony vargenes_vst')
parser$add_argument('pxc_file', metavar='pxc_file', help='peaks x cells file')
parser$add_argument('--peak_acc_crit', metavar='peak_acc_crit', default=0.005, type='double', help='peak accessible in at least X cells (if integer) or X proportion of cells (if decimal); default 0.005')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file for sample')
parser$add_argument('sample_col', metavar='sample_col', help='sample column in metadata file')
parser$add_argument('--NObinarize', action='store_true', help='if given, do NOT binarize matrix')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('--numPeaks', metavar='numPeaks', default=500, type='double', help='number of variables peaks to ask for (to start); default 500')
parser$add_argument('--step_crit', metavar='step_crit', default=50, type='double', help='step criterion; if integer - number of peaks to decrease; if decimal - percentage of peaks to decrease; default 50')
parser$add_argument('--end_crit', metavar='end_crit', default=0.8, type='double', help='end criterion; variable peaks less than X proportion of cells; default 0.8')
parser$add_argument('--recurSize', metavar='recurSize', default=200000, type='integer', help='number of columns for scale to recurse; default 200000') #very large to avoid doing this normally
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

pxc_file <- args$pxc_file
peak_acc_crit <- args$peak_acc_crit
meta_file <- args$meta_file
sample_col <- args$sample_col
NObinarize <- args$NObinarize
seed <- args$seed
numPeaks <- args$numPeaks
step_crit <- args$step_crit
end_crit <- args$end_crit
recurSize <- args$recurSize
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(pxc_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_varFeat_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)


print_time("Load matrix")
pxc <- readRDS(pxc_file)
dim(pxc)
pxc[1:3,1:3]
ll <- max(pxc)
if(ll!=1 & !NObinarize){
    cat(paste('Max feature counts',ll,'Binarizing.\n'))
    pxc <- makeBinary(pxc)
}
nCells_tot <- ncol(pxc)


print_time("Feature Accessibility Cutoff")
#this is looking at cells, not fragments!
if(NObinarize){
    pxc_rs <- rowSums(makeBinary(pxc))
} else {
    pxc_rs <- rowSums(pxc)
}
length(pxc_rs)
head(pxc_rs)
summary(pxc_rs)

cat(paste('\nFeatures accessibility cutoff:',peak_acc_crit,'\n'))
if(peak_acc_crit<1) nCells_cutoff <- floor(nCells_tot * peak_acc_crit) else nCells_cutoff <- peak_acc_crit
cat(paste('Features must be accessible in at least',nCells_cutoff,'cells\n'))

peaks_keeping <- names(pxc_rs[pxc_rs >= nCells_cutoff])
cat(paste('Features kept:',length(peaks_keeping),'\n\n'))
saveRDS(peaks_keeping,paste(sep='',outPrefix,'_accFeat-',nCells_cutoff,'.rds'))

pxc <- pxc[peaks_keeping,]
dim(pxc)
pxc[1:3,1:3]
nPeaks_tot <- nrow(pxc)


print_time("Load meta")
meta <- readRDS(meta_file)
dim(meta)
head(meta)
if(!(sample_col %in% colnames(meta))) stop('sample column not in meta file')
if(!all(colnames(pxc) %in% rownames(meta))) stop('Not all mat cells in meta cells')
meta <- meta[colnames(pxc),]
if(!identical(colnames(pxc),rownames(meta))) stop('cellnames not the same between mat and meta files')


print_time("Normalize")
set.seed(seed)
pxc_norm <- TF.IDF(pxc)
pxc_norm <- LogNormalize(pxc_norm)
saveRDS(pxc_norm,paste(sep='',outPrefix,'_norm.rds'))


print_time("Variable")
nPeaks_end <- floor(nCells_tot * end_crit)
if(step_crit<1) nPeaks_step <- floor(nPeaks_tot * step_crit)*-1 else nPeaks_step <- step_crit*-1
cat(paste("Total number of (subsetted) features:",nPeaks_tot,"\n"))
cat(paste("Ending feature criteria:",nPeaks_end,"\n"))
cat(paste("Feature Step:",nPeaks_step,"\n"))
cat(paste("Initially asking for:",numPeaks,"\n"))
for(nPeaks_ask in seq(numPeaks,1,nPeaks_step)){
    
    set.seed(seed)
    var_peaks = vargenes_vst(pxc_norm, groups = meta[,sample_col], topn = nPeaks_ask)
    saveRDS(var_peaks,paste(sep='',outPrefix,'_varFeat.rds'))
    nPeaks_var <- length(var_peaks)
    print_time(paste('Asked for',nPeaks_ask,'Got',nPeaks_var,'\n'))
    
    if(nPeaks_var < nPeaks_end) break
}


print_time("Subset, scale, and transpose") 
#want features scaled to mean 0 and var 1
cxp_scaled <- scale_byHalf(t(pxc_norm[var_peaks,]),num_cols_max=recurSize)
dim(cxp_scaled)
cxp_scaled[1:3,1:3]
mean(cxp_scaled[,1])
var(cxp_scaled[,1])
identical(var_peaks,colnames(cxp_scaled))
identical(rownames(cxp_scaled),rownames(meta))
saveRDS(cxp_scaled,paste(sep='',outPrefix,'_norm_scaled_transpose.rds'))


print_time("Done.")

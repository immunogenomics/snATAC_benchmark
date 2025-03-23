print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(magrittr)
    library(data.table)
    library(tibble)
    library(dplyr)
    library(argparse)
})


NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
    A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    A@x <- scaling_factor * A@x
    if (do_ftt) {
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else {
        A@x <- log(1 + A@x)
    }
    return(A)
}

FindVariableGenesBatch <- function(exprs_mat, meta_df, cell_column, sample_column, genes_exclude=NULL, ngenes_use=1e3, expr_min=0.1) {
    if (!is.null(genes_exclude)) {
        genes_use <- setdiff(row.names(exprs_mat), genes_exclude)
    }
    x_res <- split(meta_df[,cell_column], meta_df[,sample_column]) %>% lapply(function(x) {
        FindVariableGenesSeurat(exprs_mat[genes_use, x]) %>% 
            subset(gene.mean >= expr_min) %>% 
            tibble::rownames_to_column("gene") %>% 
            dplyr::arrange(-gene.dispersion) %>%
            head(ngenes_use)
    })
    data.table(Reduce(rbind, x_res))[, .N, by = gene][order(-N)]    
}

FindVariableGenesSeurat <- function (data, x.low.cutoff = 0.1, x.high.cutoff = 8,
                                     y.cutoff = 1, y.high.cutoff = Inf, num.bin = 0,
                                     binning.method = "equal_width", sort.results = TRUE,
                                     display.progress = TRUE, ...)
{
    genes.use <- rownames(data)
    if (class(data) != "dgCMatrix") {
        data <- as(as.matrix(data), "dgCMatrix")
    }
    ## (1) get means and variances
    gene.mean <- FastExpMean(data, display.progress)
    names(gene.mean) <- genes.use
    gene.dispersion <- FastLogVMR(data, display.progress)
    names(gene.dispersion) <- genes.use

    gene.dispersion[is.na(x = gene.dispersion)] <- 0
    gene.mean[is.na(x = gene.mean)] <- 0

    mv.df <- data.frame(gene.mean, gene.dispersion)
    rownames(mv.df) <- rownames(data)

    ## (OPTIONAL) do the binning correction
    if (num.bin > 0) {
      if (binning.method == "equal_width") {
          data_x_bin <- cut(x = gene.mean, breaks = num.bin)
      }
      else if (binning.method == "equal_frequency") {
          data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean >
              0], probs = seq(0, 1, length.out = num.bin))))
     }
      else {
          stop(paste0("Invalid selection: '", binning.method,
              "' for 'binning.method'."))
      }
      names(x = data_x_bin) <- names(x = gene.mean)
      mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = mean)
      sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = sd)
      gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
      gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0

      mv.df$gene.dispersion.scaled <- gene.dispersion.scaled
    }

    return(mv.df)
}
environment(FindVariableGenesSeurat) <- asNamespace("Seurat")

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


parser <- ArgumentParser(description='norm, variable genes, scale genes x cells file using modified Seurat methods')
parser$add_argument('gxc_file', metavar='gxc_file', help='genes x cells file')
parser$add_argument('--is_norm', action='store_true', help='if given, mat file is already normalized')
parser$add_argument('--norm_scale_factor', metavar='norm_scale_factor', default=1e4, type='double', help='scaling factor for normalization; default=10000')
parser$add_argument('--exp_min', metavar='exp_min', default=0.1, type='double', help='expression minimum; default 0.1')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file for sample')
parser$add_argument('sample_col', metavar='sample_col', help='sample column in metadata file')
parser$add_argument('--cell_col', metavar='cell_col', help='cell column in metadata file; if not given, rownames assumed')
parser$add_argument('--seed', metavar='seed', default=0, type='integer', help='randomization seed; default 0')
parser$add_argument('--numGenes', metavar='numGenes', default=500, type='double', help='number of variable genes; default 500')
parser$add_argument('--recurSize', metavar='recurSize', default=200000, type='integer', help='number of columns for scale to recurse; default 200000') #very large to avoid doing this normally
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

gxc_file <- args$gxc_file
is_norm <- args$is_norm
norm_scale_factor <- args$norm_scale_factor
exp_min <- args$exp_min
meta_file <- args$meta_file
sample_col <- args$sample_col
cell_col <- args$cell_col
seed <- args$seed
numGenes <- args$numGenes
recurSize <- args$recurSize
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(gxc_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_varGenes_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)


print_time("Load files")
meta <- readRDS(meta_file)
if(!(sample_col %in% colnames(meta))) stop('Sample column must be in metadata file')
if(is.null(cell_col)){
    cell_col <- 'cell'
    while(cell_col %in% colnames(meta)){
        cell_col <- paste(sep='',cell_col,1)
    }
    meta[,cell_col] <- rownames(meta)
}
if(!(cell_col %in% colnames(meta))) stop('Cell column must be in metadata file')

gxc <- readRDS(gxc_file)
dim(gxc)
gxc[1:3,1:3]

if(!identical(colnames(gxc),rownames(meta))){
    if(!identical(sort(colnames(gxc)),sort(rownames(meta)))){
        stop('cell names between gxc and meta must agree.')
    } else {
        meta <- meta[colnames(gxc),]
    }
}


print_time("Normalize matrix")
if(!is_norm){
    gxc_norm <- gxc %>% NormalizeDataSeurat(scaling_factor=norm_scale_factor)
    saveRDS(gxc_norm,paste(sep='',outPrefix,'_norm.rds'))
} else {
    gxc_norm <- gxc[1:nrow(gxc),1:ncol(gxc)]
}


print_time("Gene QC for vargenes")
genes_toExclude <- grep("^MT-|^RPL|^RPS|MALAT1|^MIR", row.names(gxc_norm), value = TRUE)
cat(paste("Removing",length(genes_toExclude),"MT, ribosomal, and miRNA genes.\n"))


print_time("Find variable genes")
vargenes_df <- FindVariableGenesBatch(gxc_norm, meta, cell_col, sample_col,
                                      genes_exclude=genes_toExclude, ngenes_use=numGenes, expr_min=exp_min)
var_genes <- vargenes_df$gene


print_time("Subset, scale, and transpose")
#want features scaled to mean 0 and var 1
cxg_scaled <- scale_byHalf(t(gxc_norm[var_genes, ]),num_cols_max=recurSize)
dim(cxg_scaled)
cxg_scaled[1:3,1:3]
mean(cxg_scaled[,1])
var(cxg_scaled[,1])
identical(var_genes,colnames(cxg_scaled))
identical(rownames(cxg_scaled),rownames(meta))
saveRDS(cxg_scaled, paste(sep='',outPrefix,'_norm_scaled_transpose.rds'))


print_time('Done.')

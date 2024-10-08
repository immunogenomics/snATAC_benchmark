print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(stringr)
    library(lisi)
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

regular_lisi_function <- function(mat,met,cols,addCellCol=FALSE){
    if(!all(cols %in% colnames(met))) stop('not all cols in meta')
    if(!identical(rownames(met),rownames(mat))) stop('not same cells/order')
    
    ret_lisi_df <- rbind(lisi::compute_lisi(mat, meta, cols))
    print(colMeans(ret_lisi_df))
    
    if(addCellCol){
        ret_lisi_df$cell <- rownames(ret_lisi_df)
    }
    
    return(ret_lisi_df)
}

skip_lisi_function <- function(mat,met,col,skipVal,addCellCol=FALSE){
    if(length(col)!=1) stop('it seems unlikely that the same cells will have the same missing values, so do this function one at a time.')
    if(!(col %in% colnames(met))) stop('not all col in meta')
    if(!any(skipVal %in% met[,col])) stop(paste('no',skipVal,'in column',col,', so use regular_lisi_function'))
    if(!identical(rownames(met),rownames(mat))) stop('not same cells/order')
    
    wSkip_cells <- rownames(met[which(!(met[,col] %in% skipVal)),])
    ret_lisi_df <- lisi::compute_lisi(mat[wSkip_cells,], met[wSkip_cells,], col)
    print(colMeans(ret_lisi_df))
    
    if(addCellCol){
        ret_lisi_df$cell <- rownames(ret_lisi_df)
    }
    
    return(ret_lisi_df)
}

skip_lisi_iter_function <- function(mat,met,cols,skipVal){
    if(!all(cols %in% colnames(met))) stop('not all cols in meta')
    if(!identical(rownames(met),rownames(mat))) stop('not same cells/order')
    
    if(length(cols)>0){
        ret_lisi_df <- skip_lisi_function(mat,met,cols[1],skipVal,addCellCol=TRUE)
    }
    if(length(cols)>1){
        for(cc in cols[2:length(cols)]){
            ll <- skip_lisi_function(mat,met,cc,skipVal,addCellCol=TRUE)
        }
        ret_lisi_df <- merge(ret_lisi_df,ll,all=TRUE)
        rownames(ret_lisi_df) <- ret_lisi_df$cell
    }
    
    return(ret_lisi_df)
}



parser <- ArgumentParser(description='snATAC-seq benchmark LISI-based metrics')
parser$add_argument('snATAC_embed_file', metavar='snATAC_embed_file', help='snATAC-seq embedding file')
parser$add_argument('meta_file', metavar='meta_file', help='meta file - all embed rownames should be in this file')
parser$add_argument('lisi_cols', metavar='lisi_cols', help='columns for lisi, separated by commas')
parser$add_argument('--skip_val_lisi', metavar='skip_val_lisi', help='Skip these values when calculating LISI, separated by commas')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

snATAC_embed_file <- args$snATAC_embed_file
meta_file <- args$meta_file
lisi_cols <- args$lisi_cols
skip_val_lisi <- args$skip_val_lisi
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(snATAC_embed_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time('Load files')
embMat <- readRDS(snATAC_embed_file)
meta <- readRDS(meta_file)

embMat <- as.matrix(embMat)

if(!any(rownames(embMat) %in% rownames(meta))){
    cat('Embedding matrix cell names did not match meta matrix cell names; trying ArchR style.\n')
    split1 <- str_split_fixed(rownames(embMat),'#',2)
    split2 <- str_split_fixed(split1[,2],'-',2)
    rownames(embMat) <- paste(sep='',split1[,1],'_',split2[,1])
}
if(!all(rownames(embMat) %in% rownames(meta))) stop('All embedding rownames need to be in meta rownames.')
meta <- meta[rownames(embMat),]

lisi_col_toUse <- parse_columns(lisi_cols)
if(!all(lisi_col_toUse %in% colnames(meta))) stop('Not all lisi columns in meta.')


print_time("LISI")
if(!is.null(skip_val_lisi)){
    skip_val_lisi <- parse_columns(skip_val_lisi)
    
    cols_logic <- apply(meta[,lisi_col_toUse],2,FUN=function(x){any(grepl(paste(skip_val_lisi,collapse="|"),x))})
    cols_wSkip <- names(cols_logic[cols_logic==TRUE])
    cols_noSkip <- names(cols_logic[cols_logic==FALSE])
    
    if(length(cols_noSkip)>0 & length(cols_wSkip)>0){
        #do both
        lisi_noSkip_df <- regular_lisi_function(embMat,meta,cols_noSkip,addCellCol=TRUE)
        lisi_wSkip_df <- skip_lisi_iter_function(embMat,meta,cols_wSkip,skip_val_lisi)
        lisi_df <- merge(lisi_noSkip_df,lisi_wSkip_df,all.x=TRUE)
        rownames(lisi_df) <- lisi_df$cell
        lisi_df <- lisi_df[,which(colnames(lisi_df)!='cell')]
    } else if(length(cols_noSkip)>0 & length(cols_wSkip)==0){
        #do noskip
        lisi_df <- regular_lisi_function(embMat,meta,cols_noSkip,addCellCol=FALSE)
    } else if(length(cols_noSkip)==0 & length(cols_wSkip)>0){
        #do skip
        lisi_df <- skip_lisi_iter_function(embMat,meta,cols_wSkip,skip_val_lisi)
    } else {
        stop('no LISI columns')
    }
} else {
    lisi_df <- regular_lisi_function(embMat,meta,lisi_col_toUse,addCellCol=FALSE)
}
colnames(lisi_df) <- paste(sep='','lisi_',colnames(lisi_df))
saveRDS(lisi_df,paste(sep='',outPrefix,'_LISImetrics.rds'))


print_time('Done.')


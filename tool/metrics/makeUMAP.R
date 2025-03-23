print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(umap)
    library(ggplot2)
    library(ggrastr)
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


plot_dot <- function(toPlot, plotX, plotY, plotTit, plotPrefix,
                     plotHeight=6, plotWidth=6, plotMeta=FALSE, plotMetaCols=NA,
                     plotColors=NA, numeric_color_low='palegreen1', numeric_color_high='palegreen4'){

    if(!all(c(plotX, plotY) %in% colnames(toPlot))) stop('columns not in df')

    g <- ggplot(toPlot,aes(x=!!sym(plotX),y=!!sym(plotY))) + rasterise(geom_point(size=1,alpha=0.5),dpi=300) +
            theme_classic(base_size = 20) + ggtitle(plotTit)
    ggsave(paste(sep="",plotPrefix,".png"),plot=g,units='in',height=plotHeight,width=plotWidth)
    
    if(plotMeta){
        for(col in plotMetaCols){
            g <- ggplot(toPlot,aes(x=!!sym(plotX),y=!!sym(plotY),color=!!sym(col))) + 
                    geom_point(size=1,alpha=0.5) + theme_classic(base_size = 20) + 
                    ggtitle(plotTit)
            
            if(is.numeric(toPlot[,col])){
                g <- g + scale_color_gradient(low=numeric_color_low,high=numeric_color_high)
                
                widthAdd = floor(nchar(col)/7)
            } else {
                if(all(unique(toPlot[,col]) %in% names(plotColors))){
                    g <- g + scale_color_manual(values=plotColors) + 
                            guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))
                }
                
                widthAdd = floor(max(nchar(unique(toPlot[,col])),nchar(col))/7)
                if(length(unique(toPlot[,col]))>17) widthAdd=widthAdd*2
            }
            ggsave(paste(sep="",plotPrefix,"_col-",col,".png"),
                   plot=g,units='in',height=plotHeight,width=plotWidth+widthAdd)
        }
    }
}


parser <- ArgumentParser(description='UMAP')
parser$add_argument('emb_file', metavar='emb_file', help='cell x dim file')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('--meta_cols', metavar='meta_cols', help='columns to plot in meta file; comma delimited list; if not given, only 1 UMAP plotted.')
parser$add_argument('--LISI_file', metavar='LISI_file', help='if given, add to meta and subset by meta_cols')
parser$add_argument('--NN_file', metavar='NN_file', help='if given, add to meta and subset by meta_cols')
parser$add_argument('plotTitle', metavar='plotTitle', help='title for plots')
parser$add_argument('--color_file', metavar='color_file', help='file for colors')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='integer', help='randomization seed; default 1234567890')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

emb_file <- args$emb_file
meta_file <- args$meta_file
meta_cols <- args$meta_cols
LISI_file <- args$LISI_file
NN_file <- args$NN_file
plotTitle <- args$plotTitle
color_file <- args$color_file
seed <- args$seed
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(emb_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_UMAP_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Overall seed")
set.seed(seed)


print_time('Load embedding file')
embMat <- readRDS(emb_file)
dim(embMat)
embMat[1:3,1:3]


print_time('Load meta file')
meta <- readRDS(meta_file)

if(!is.null(LISI_file)){
    if(!file.exists(LISI_file)){
        cat("SKIPPING: LISI file does not exist.\n")
    } else {
        cat("Adding LISI file to meta.\n")
        LISI_df <- readRDS(LISI_file)
        if(!all(rownames(LISI_df) %in% rownames(meta))) stop('Cell names in LISI do not match metadata')
        meta <- cbind(meta,LISI_df[rownames(meta),])
    }
}

if(!is.null(NN_file)){
    if(!file.exists(NN_file)){
        cat("SKIPPING: NN file does not exist.\n")
    } else {
        cat("Adding NN file to meta.\n")
        NN_df <- readRDS(NN_file)
        if(!all(rownames(NN_df) %in% rownames(meta))) stop('Cell names in NN do not match metadata')
        meta <- cbind(meta,NN_df[rownames(meta),])
    }
}


print_time("Verify cell names")
if(!any(rownames(embMat) %in% rownames(meta))){
    #fix embedding matrix cell names under the assumption that they are ArchR style, then verify
    split1 <- str_split_fixed(rownames(embMat),'#',2)
    split2 <- str_split_fixed(split1[,2],'-',2)
    cellNames <- paste(sep='',split1[,1],'_',split2[,1])
    print(head(c(cellNames,rownames(embMat))))
    
    if(!any(cellNames %in% rownames(meta))) stop('Cell names in embedding matrix do not match metadata or ArchR style')
    rownames(embMat) <- cellNames
}


print_time("Meta columns")
meta_col_toUse <- parse_columns(meta_cols)

use_meta <- FALSE
if(all(rownames(embMat) %in% rownames(meta))){
    meta <- meta[rownames(embMat),]

    if(any(meta_col_toUse %in% colnames(meta))){
        meta_col_toUse <- meta_col_toUse[which(meta_col_toUse %in% colnames(meta))]
        meta <- meta[,meta_col_toUse,drop=FALSE]
        use_meta <- TRUE
    }
}


print_time("Colors")
color_vec <- NA
if(use_meta & !is.null(color_file)){
    if(file.exists(color_file)){
        print_time('colors file')
        color_vec <- readRDS(color_file)
    }
}


print_time('UMAP')
set.seed(seed)
umap_res <- umap(embMat, n_neighbors = 30, metric = "cosine", min_dist = .3) 
saveRDS(umap_res,paste(sep='',outPrefix,'_UMAP_res.rds'))
umap_df <- data.frame("UMAP1"=umap_res$layout[, 1],"UMAP2"=umap_res$layout[, 2])
saveRDS(umap_df,paste(sep='',outPrefix,'_UMAP.rds'))


print_time('Full DF')
full_df <- cbind(meta,embMat,umap_df)
set.seed(seed)
full_df <- full_df[sample(nrow(full_df),nrow(full_df)),]


print_time('Plot initial dimensions')
plot_dot(full_df,colnames(embMat)[1],colnames(embMat)[2],plotTitle, paste(sep="",outPrefix,'_dim1-dim2'), 
         plotMeta=use_meta, plotMetaCols=meta_col_toUse, plotColors=color_vec)


print_time("Plot UMAPs")
plot_dot(full_df,'UMAP1','UMAP2',paste(plotTitle,ncol(embMat),'dim'), paste(sep="",outPrefix,'_UMAP'), 
         plotMeta=use_meta, plotMetaCols=meta_col_toUse, plotColors=color_vec)


print_time('Done.')


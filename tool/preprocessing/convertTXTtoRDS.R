print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(stringr)
    library(plyr)
    library(argparse)
})


parser <- ArgumentParser(description='convert matrix from txt to rds')
parser$add_argument('txt_file', metavar='txt_file', help='matrix txt file with row and column names')
parser$add_argument('--addHyp', action='store_true', help='if given, add hyphen back')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file to verify cellnames')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

txt_file <- args$txt_file
addHyp <- args$addHyp
meta_file <- args$meta_file
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(txt_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_convRDS_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load files")
meta <- readRDS(meta_file)
nrow(meta)
mat <- read.table(txt_file,header=TRUE,row.names=1,stringsAsFactors=FALSE)
nrow(mat)
if(addHyp){
    cellName_conv_df <- data.frame('wHyp'=rownames(meta),'nHyp'=gsub('-','',rownames(meta)),stringsAsFactors=FALSE)
    new_mat_cellName <- mapvalues(rownames(mat),from=cellName_conv_df$nHyp,to=cellName_conv_df$wHyp)
    print(head(cbind(rownames(mat),new_mat_cellName)))
    rownames(mat) <- new_mat_cellName
}
dim(mat)
mat[1:3,1:3]


print_time("Verify cell order")
if(!identical(sort(rownames(meta)),sort(rownames(mat)))) stop('cells not identical')
mat <- mat[rownames(meta),]


print_time("Save")
saveRDS(mat,paste(sep='',outPrefix,'.rds'))


print_time('Done.')


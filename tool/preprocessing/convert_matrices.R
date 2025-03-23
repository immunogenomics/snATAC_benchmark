print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(stringr)
    library(Biostrings)
    library(dplyr)
    library(gtools)
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

#COPIED from Signac::BinarizeCounts
makeBinary <- function(object) {
    if (inherits(x = object, what = "CsparseMatrix")) {
        slot(object = object, name = "x") <- rep.int(x = 1, times = length(x = slot(object = object, name = "x")))
    } else {
        object[object > 1] <- 1
    }
    return(object)
}


parser <- ArgumentParser(description='convert RDS Matrix to 10x MM matrix with cell and feature files')
parser$add_argument('mat_file', metavar='mat_file', help='feature x cells matrix file')
parser$add_argument('--rmHyp', action='store_true', help='if given, remove hyphen')
parser$add_argument('--binarize', action='store_true', help='if given, binarize matrix')
parser$add_argument('--varFeat_file', metavar='varFeat_file', help='variable features file: RDS with vector of chr:start-stop')
parser$add_argument('--mat_transpose', action='store_true', help='if given, save transposed matrix, i.e., cells x feature')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file')
parser$add_argument('--writeCells', action='store_true', help='if given, write cell names to file')
parser$add_argument('--conv_meta_txt', action='store_true', help='if given, convert full metadata file from RDS to TXT')
parser$add_argument('--meta_cols', metavar='meta_cols', help='columns to include in barcodes.tsv; comma delimited list; if not given, no columns.')
parser$add_argument('--feat_form', metavar='feat_form', default='tsv', choices=c('tsv','bed','hg19','hg38'), help='Output feature file format; tsv is chr:start-stop; bed is chr<tab>start<tab>stop; hg19 and hg38 are genome versions for FASTA files; default tsv')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
args <- parser$parse_args()

mat_file <- args$mat_file
rmHyp <- args$rmHyp
binarize <- args$binarize
varFeat_file <- args$varFeat_file
mat_transpose <- args$mat_transpose
meta_file <- args$meta_file
writeCells <- args$writeCells
conv_meta_txt <- args$conv_meta_txt
meta_cols <- args$meta_cols
feat_form <- args$feat_form
outDir <- args$outDir

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Load files")
mat <- readRDS(mat_file)
dim(mat)
mat[1:3,1:3]
ll <- max(mat)
if(ll!=1 & binarize){
    mat <- makeBinary(mat)
    cat(paste('Binarized. New max:',max(mat),'\n'))
}

if(!is.null(varFeat_file)){
    varFeat <- readRDS(varFeat_file)
    varFeat <- sort(varFeat)
    length(varFeat)
    head(varFeat,n=3)
    if(!all(varFeat %in% rownames(mat))) stop('not all varFeat in matrix rownames')
    
    mat <- mat[varFeat,]
}

meta <- readRDS(meta_file)
nrow(meta)
head(rownames(meta),n=3)
if(!all(rownames(meta) %in% colnames(mat))) stop('not all cells in matrix colnames')
mat <- mat[,rownames(meta)]

if(rmHyp){
    colnames(mat) <- str_replace(colnames(mat),'-','')
    rownames(meta) <- str_replace(rownames(meta),'-','')
}

print_time("Write Files")
dim(mat)
mat[1:3,1:3]

if(feat_form=='bed'){
    split <- str_split_fixed(rownames(mat),':',2)
    split2 <- str_split_fixed(split[,2],'-',2)
    peaks <- data.frame('chr'=split[,1],'start'=split2[,1],'stop'=split2[,2])
    write.table(peaks,file=paste(sep='',outDir,'peaks.bed'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if(feat_form=='tsv'){
    write(rownames(mat),ncolumns=1,file=paste(sep='',outDir,'features.tsv'))
} else if(feat_form %in% c('hg19','hg38')){
    print_time('Make Feature GRanges')
    split <- str_split_fixed(rownames(mat),':',2)
    split2 <- str_split_fixed(split[,2],'-',2)
    gr <- GRanges(seqnames = split[,1],
                          ranges = IRanges(start = as.numeric(split2[,1]), end = as.numeric(split2[,2]),
                                           names = rownames(mat)))
    print(gr)
    if(!identical(rownames(mat),names(gr))) stop('features not identical')
    
    
    print_time("Get FASTA")
    if(feat_form=='hg38'){
        gen <- BSgenome.Hsapiens.UCSC.hg38
    } else if(feat_form=='hg19'){
        gen <- BSgenome.Hsapiens.UCSC.hg19
    } else {
        stop('Should not have gotten here.')
    }
    getSeq(gr, x=gen) %>% writeXStringSet(filepath=paste(sep='',outDir,'features.fa'))

} else {
    stop('Should not have gotten here.')
}

if(writeCells){
    write(colnames(mat),ncolumns=1,file=paste(sep='',outDir,'cells.txt'))
}
    
if(conv_meta_txt){
    write.table(meta,file=paste(sep='',outDir,'meta.txt'),sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
}

if(!is.null(meta_cols)){
    meta_cols <- parse_columns(meta_cols)
    if(!any(meta_cols %in% colnames(meta))){
        cat("None of the meta columns given are in this meta.\n")
        meta_cols <- NULL
        
    } else {
        meta_cols <- meta_cols[which(meta_cols %in% colnames(meta))]
        meta <- meta[,meta_cols]
        print(head(meta,n=3))

        write.table(meta,file=paste(sep='',outDir,'barcodes.tsv'),sep='\t',row.names=TRUE,col.names=FALSE,quote=FALSE)
        write(colnames(meta),ncolumns=1,file=paste(sep='',outDir,'meta_colnames.tsv'))
    }
}
if(is.null(meta_cols)){
    write(colnames(mat),ncolumns=1,file=paste(sep='',outDir,'barcodes.tsv'))
}

if(mat_transpose){
    writeMM(t(mat),file=paste(sep='',outDir,'matrix_transpose.mtx'))
} else {
    writeMM(mat,file=paste(sep='',outDir,'matrix.mtx'))
}

print_time("Done.")

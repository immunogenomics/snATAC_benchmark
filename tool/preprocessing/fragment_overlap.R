print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(GenomicRanges)
    library(Matrix)
    library(Matrix.utils)
    library(gtools)
    library(argparse)
})

gRange_str <- function(gr) {
    return(paste(sep="",seqnames(gr),":",start(gr),"-",end(gr)))
}


parser <- ArgumentParser(description='overlap feature by cell fragments to get matrix')
parser$add_argument('feature_file', metavar='feature_file', help='feature file; no header - usually peaks or genes')
parser$add_argument('--feat_colNum', metavar='feat_colNum', type='integer', help='if given, feature file column number to use as feature name; if not given, default to col1:col2-col3')
parser$add_argument('fragment_file', metavar='fragment_file', help='fragment file; no header')
parser$add_argument('--donor_name', metavar='donor_name', default='donor', help='donor for use in cell name; default donor')
parser$add_argument('--donor_delim', metavar='donor_delim', default='_', help='delimiter between donor and cell barcode; default _')
parser$add_argument('--NOrm1', action='store_true', help='if given, do NOT remove the "-1" at the end of the cell barcode; note it will verify that all cell barcodes have it before removing any!')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

feature_file <- args$feature_file
feat_colNum <- args$feat_colNum
fragment_file <- args$fragment_file
donor_name <- args$donor_name
donor_delim <- args$donor_delim
NOrm1 <- args$NOrm1
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(feature_file,fragment_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_fragOverlap_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Staring feature GR")
peaks_df <- read.table(feature_file,sep="\t",stringsAsFactors=FALSE)
if(ncol(peaks_df)<3) stop('Feature file should have at least 3 columns: chr, start, end')
if(!is.null(feat_colNum)){ 
    if(ncol(peaks_df)<feat_colNum) stop('Too few columns for given feat_colNum.') 
    if(any(duplicated(peaks_df[,feat_colNum]))) stop('Name column in feature file cannot have duplicates.')
}
full_colnames <- c('chr','start','end')
if(ncol(peaks_df)>3) {full_colnames <- c(full_colnames,paste(sep='','col',4:ncol(peaks_df)))}
colnames(peaks_df)<-full_colnames[1:ncol(peaks_df)]
if(!is.null(feat_colNum)) {colnames(peaks_df)[feat_colNum] <- 'geneName'}
head(peaks_df)
peaks_GR <- GRanges(peaks_df)
peaks_GR


print_time("Starting fragment GR")
fragment_df <- read.table(fragment_file,sep="\t",stringsAsFactors=FALSE)
if(ncol(fragment_df)<4) stop('Fragment file should have at least 4 columns: chr, start, end, cell')
full_colnames<-c('chr','start','end','cell','dupCt')
if(ncol(fragment_df)>5) {full_colnames <- c(full_colnames,paste(sep='','col',6:ncol(fragment_df)))}
colnames(fragment_df)<-full_colnames[1:ncol(fragment_df)]
fragment_df$donor <- donor_name
fragment_df$CBD <- paste(sep=donor_delim,fragment_df$donor,fragment_df$cell)
if(!NOrm1){ 
    if(all(substrRight(fragment_df$CBD,2)=='-1')){
        fragment_df$CBD <- substr(fragment_df$CBD,1,nchar(fragment_df$CBD)-2)
    }
}
head(fragment_df)
fragment_GR <- GRanges(fragment_df)
fragment_GR
length(unique(fragment_GR$cell))


print_time("Starting initial findOverlaps for reordering")
fragment_GR_forReorder <- fragment_GR[seqnames(fragment_GR) == as.character(seqnames(fragment_GR)[1])]
peaks_GR_forReorder <- peaks_GR[seqnames(peaks_GR) == as.character(seqnames(fragment_GR)[1])]
hits_init <- findOverlaps(fragment_GR_forReorder,peaks_GR_forReorder)
peaks_GR <- c(peaks_GR[peaks_GR != peaks_GR_forReorder[subjectHits(hits_init)[1]]],peaks_GR[peaks_GR == peaks_GR_forReorder[subjectHits(hits_init)[1]]])
peaks_GR
fragment_GR <- c(fragment_GR[fragment_GR != fragment_GR_forReorder[queryHits(hits_init)[1]]],fragment_GR[fragment_GR == fragment_GR_forReorder[queryHits(hits_init)[1]]])
fragment_GR
rm(fragment_GR_forReorder,peaks_GR_forReorder)
gc(verbose=FALSE)


print_time("Starting findOverlaps")
hits <- findOverlaps(fragment_GR,peaks_GR)
hits_SM <- sparseMatrix(i=queryHits(hits),j=subjectHits(hits),dimnames = list(gRange_str(fragment_GR),gRange_str(peaks_GR)))
hits_SM[1:3,1:3]
dim(hits_SM)


print_time("Starting peaks by cells")
rownames(hits_SM)<-fragment_GR$CBD
cells_by_peaks <- aggregate.Matrix(hits_SM, row.names(hits_SM))
peaks_by_cells <- t(cells_by_peaks)
rm(cells_by_peaks)
gc(verbose=FALSE)
if(!is.null(feat_colNum)){ rownames(peaks_by_cells) <- peaks_GR$geneName }
peaks_by_cells <- peaks_by_cells[mixedsort(rownames(peaks_by_cells)),]
peaks_by_cells[1:3,1:3]
dim(peaks_by_cells)
saveRDS(peaks_by_cells,paste(sep="",outDir,prefix,"Xcells.rds"))


print_time("Starting peaks by donors")
rs <- rowSums(peaks_by_cells)
peaksXdonor <- data.frame(V1=ifelse(rs==0,0,1))
peaksXdonor$V1 <- as.integer(peaksXdonor$V1)
colnames(peaksXdonor) <- donor_name
table(peaksXdonor)
saveRDS(peaksXdonor,paste(sep="",outDir,prefix,"Xdonor.rds"))


print_time("Done.")


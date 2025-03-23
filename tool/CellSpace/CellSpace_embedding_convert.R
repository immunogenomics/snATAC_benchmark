print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(CellSpace)
    library(Matrix)
    library(argparse)
})


parser <- ArgumentParser(description='convert CellSpace cell embeddings from tsv to rds')
parser$add_argument('cellspace_tsv_file', metavar='cellspace_tsv_file', help='cellspace tsv file')
parser$add_argument('--cellspace_name', metavar='cellspace_name', default='project', help='cellspace project name; default: project')
parser$add_argument('cellspace_cells_file', metavar='cellspace_cells_file', help='cellspace cells file')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

cellspace_tsv_file <- args$cellspace_tsv_file
cellspace_name <- args$cellspace_name
cellspace_cells_file <- args$cellspace_cells_file
meta_file <- args$meta_file
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(cellspace_tsv_file,meta_file,cellspace_cells_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_convCSemb_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load cell files")
cell.idx <- readLines(cellspace_cells_file)
length(cell.idx)
meta <- readRDS(meta_file)
nrow(meta)


print_time("Verify cell order")
if(!all(cell.idx %in% rownames(meta))) stop('not all cells from cellspace in meta')
meta <- meta[cell.idx,]
if(!identical(rownames(meta),cell.idx)) stop('cells not identical')


print_time("Make Cell Space object")
cso <- CellSpace(
  project = cellspace_name,
  emb.file = cellspace_tsv_file,
  meta.data = meta
)
saveRDS(cso,paste(sep='',outPrefix,'_CS_object.rds'))


print_time("Get Cell Space Cell Embeddings")
dim(cso@cell.emb)
cso@cell.emb[1:3,1:3]
saveRDS(cso@cell.emb,paste(sep='',outPrefix,'_CS_embeddings.rds'))


print_time('Done.')


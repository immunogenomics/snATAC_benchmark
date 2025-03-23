print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(SummarizedExperiment)
    library(Matrix)
    library(stringr)
    library(stringi)
    library(argparse)
})

convert_cellNames <- function(cb,delimiter='_'){
    spl <- str_split_fixed(stri_reverse(cb),'_',2)
    ret <- paste(sep='#', stri_reverse(spl[,2]),stri_reverse(spl[,1]))
    if(all(substrRight(ret,2)!='-1')) ret <- paste(sep='',ret,'-1')
    return(ret)
}


parser <- ArgumentParser(description='Make basicGAS SummarizedExperiment for ArchR Input')
parser$add_argument('mat_file', metavar='mat_file', help='matrix file - basicGAS x cells')
parser$add_argument('--delim', metavar='delim', default='_', help='cell name delimiter between sample and cell barcode; default _')
parser$add_argument('meta_file', metavar='meta_file', help='meta file - cells x metadata')
parser$add_argument('gr_file', metavar='gr_file', help='GRanges file - basicGAS')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path; cells sample#CB-1')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

mat_file <- args$mat_file
delim <- args$delim
meta_file <- args$meta_file
gr_file <- args$gr_file
project_dir <- args$project_dir
wk_dir <- args$wk_dir
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(mat_file,meta_file,gr_file,project_dir))) stop("Input(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(outDir,1)!='/') project_dir<-paste(sep="",outDir,'/')
if(!file.exists(outDir)) dir.create(outDir)

outPrefix <- paste(sep='',outDir,prefix)
cat(paste("Using:",outPrefix,'\n'))

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Set wd")
if(length(wk_dir)==0){
    wk_dir <- paste(sep='',dirname(project_dir),'/')
} else {
    if(substrRight(wk_dir,1)!='/') wk_dir<-paste(sep="",wk_dir,'/')
    if(!file.exists(wk_dir)) wk_dir <- paste(sep='',dirname(project_dir),'/')
}
cat(paste('New working dir:',wk_dir,'\n'))
setwd(wk_dir)


print_time('Load files')
mat <- readRDS(mat_file)
dim(mat)
mat[1:3,1:3]
meta <- readRDS(meta_file)
dim(meta)
meta[1:3,1:3]
genes_GR <- readRDS(gr_file)
genes_GR

if(!identical(colnames(mat),rownames(meta))){
    if(!identical(sort(colnames(mat)),sort(rownames(meta)))){
        stop('Cell names between mat and meta must agree.')
    } else {
        meta <- meta[colnames(mat),]
    }
}

if(!all(rownames(mat) %in% names(genes_GR))) stop('Not all basicGAS matrix genes in gene GRanges')
genes_GR_subset <- genes_GR[rownames(mat),]
genes_GR_subset
if(!identical(rownames(mat),names(genes_GR_subset))) stop('basicGAS matrix genes and gene GRanges not identical')


print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)
ArchR_cells <- getCellNames(proj)
length(ArchR_cells)
head(ArchR_cells)

cellNames <- convert_cellNames(colnames(mat),delimiter=delim)
if(!all(ArchR_cells %in% cellNames)) stop('Not all ArchR cells in mat cells')

rownames(meta) <- cellNames
colnames(mat) <- cellNames
meta <- meta[ArchR_cells,]
mat <- mat[,ArchR_cells]

dim(mat)
mat[1:3,1:3]
dim(meta)
meta[1:3,1:3]


print_time('Make SumExp')
if(!identical(rownames(meta),colnames(mat))) stop('meta and mat cell names not identical')
if(!identical(names(genes_GR_subset),rownames(mat))) stop('genes GR and mat genes not identical')

gxc_se <- SummarizedExperiment(list('counts'=mat),rowRanges=genes_GR_subset,colData=meta)
gxc_se


print_time("Saving")
saveRDS(gxc_se,paste(sep='',outPrefix,'_sumExp.rds'))


print_time("Done.")


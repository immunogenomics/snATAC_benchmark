print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(argparse)
})


parser <- ArgumentParser(description='Add IterativeLSI to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('--matrix_name', metavar='matrix_name', default='TileMatrix', help='Matrix to use for iterativeLSI; default: TileMatrix')
parser$add_argument('--GAS_opt', action='store_true', help='if given, use gene defaults for addIterativeLSI, binarize=F, firstSelection=var, not depthCol since GAS still based off scATAC')
parser$add_argument('--num_varGene', metavar='num_varGene', default=2000, type='integer', help='Number of variable genes to use (only used if GAS_opt given); default 2000')
parser$add_argument('--max_dim', metavar='max_dim', default=30, type='integer', help='max dimensions - note that ArchR removes dim correlated to fragment count; default: 30')
parser$add_argument('--binarize', action='store_true', help='if given, binarize matrix before itLSI - NOTE: if GAS_opt is given, binarize is always FALSE')
parser$add_argument('--LSI_name', metavar='LSI_name', default='IterativeLSI', help='name of IterativeLSI slot; default: IterativeLSI')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--extract_addToFileName', metavar='extract_addToFileName', help='if given, add to the extract filename')
parser$add_argument('--extract_mat', action='store_true', help='if given, extract matrix')
parser$add_argument('--mat_binarize', action='store_true', help='if given, extract matrix in binary form; enforced for TileMatrix')
parser$add_argument('--saveIter', action='store_true', help='if given, save itLSI iterations')
parser$add_argument('--thread_num', metavar='thread_num', default=8, type='integer', help='number of threads to use; default 8')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
parser$add_argument('--force', action='store_true', help='if given, force overwrite')
args <- parser$parse_args()

project_dir <- args$project_dir
matrix_name <- args$matrix_name
GAS_opt <- args$GAS_opt
num_varGene <- args$num_varGene
max_dim <- args$max_dim
binarize <- args$binarize
LSI_name <- args$LSI_name
extract_dir <- args$extract_dir
extract_addToFileName <- args$extract_addToFileName
extract_mat <- args$extract_mat
mat_binarize <- args$mat_binarize
saveIter <- args$saveIter
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed
force <- args$force

print_time("Argument Checking")
if(!all(file.exists(project_dir))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(extract_dir,1)!='/') project_dir<-paste(sep="",extract_dir,'/')
if(!file.exists(extract_dir)) dir.create(extract_dir)

if(!is.null(extract_addToFileName)){
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',matrix_name,'_',extract_addToFileName)
} else {
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',matrix_name)
}
if(GAS_opt){
    outPrefix <- paste(sep='',outPrefix,'_GASopt')
}
cat(paste("Using:",outPrefix,'\n'))

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_itLSI_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
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

print_time('Set seed and threads')
set.seed(seed)
addArchRThreads(threads = thread_num)

print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)
if(!(matrix_name %in% getAvailableMatrices(proj))) stop('A valid matrix is required for iterativeLSI')
if(matrix_name=='TileMatrix') {mat_binarize <- TRUE}
if(extract_mat){
    saveRDS(getMatrixFromProject(ArchRProj = proj, useMatrix=matrix_name, binarize = mat_binarize),
            paste(sep='',outPrefix,'.rds'))
}

print_time("Iterative LSI")
if(GAS_opt){
    proj <- addIterativeLSI(ArchRProj = proj, useMatrix = matrix_name, name = LSI_name, seed = seed, 
                            dimsToUse = 1:max_dim, force = force, saveIterations = saveIter,
                            binarize = FALSE, firstSelection = "var", varFeatures = num_varGene)
} else {
    proj <- addIterativeLSI(ArchRProj = proj, useMatrix = matrix_name, name = LSI_name, seed = seed, 
                            dimsToUse = 1:max_dim, binarize = binarize, force = force, saveIterations = saveIter)
}
saveRDS(getReducedDims(ArchRProj = proj, reducedDims = LSI_name),
        paste(sep='',outPrefix,'_IterativeLSI.rds'))

if(LSI_name %in% names(proj@reducedDims)){
    if('LSIFeatures' %in% names(proj@reducedDims[[LSI_name]])){
        varTile_df <- as.data.frame(proj@reducedDims[[LSI_name]]$LSIFeatures,stringsAsFactors=FALSE)
        options(scipen=15)
        write.table(varTile_df,file=paste(sep='',outPrefix,'_IterativeLSI_varFeat.txt'),
                    row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
    }
}

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")


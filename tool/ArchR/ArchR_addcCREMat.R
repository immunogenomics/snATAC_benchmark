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


parser <- ArgumentParser(description='Add ATAC feature x cells to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('feature_file', metavar='feature_file', help='ATAC feature file')
parser$add_argument('matrix_name', metavar='matrix_name', help='name for matrix slot in ArchR project')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--extract_addToFileName', metavar='extract_addToFileName', help='if given, add to the extract filename - matrix_name is already appended!')
parser$add_argument('--ceiling_num', metavar='ceiling_num', default=4, type='integer', help='max number of reads in matrix; ArchR default 4 (for peaks, but features default is 10^9)')
parser$add_argument('--binarize', action='store_true', help='if given, binarize matrix')
parser$add_argument('--thread_num', metavar='thread_num', default=8, type='integer', help='number of threads to use; default 8')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
parser$add_argument('--noNameCheck', action='store_true', help='If given, do not exit if matrix name already used') 
args <- parser$parse_args()

project_dir <- args$project_dir
feature_file <- args$feature_file
matrix_name <- args$matrix_name
extract_dir <- args$extract_dir
extract_addToFileName <- args$extract_addToFileName
ceiling_num <- args$ceiling_num
binarize <- args$binarize
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed
noNameCheck <- args$noNameCheck

print_time("Argument Checking")
if(!all(file.exists(project_dir,feature_file))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(extract_dir,1)!='/') project_dir<-paste(sep="",extract_dir,'/')
if(!file.exists(extract_dir)) dir.create(extract_dir)

if(!is.null(extract_addToFileName)){
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',matrix_name,'_',extract_addToFileName)
} else {
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',matrix_name)
}
cat(paste("Using:",outPrefix,'\n'))

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_addcCREMat_args.txt")
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

print_time("Load Features")
features_GR <- readRDS(feature_file)

print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)
if(!noNameCheck & (matrix_name %in% getAvailableMatrices(proj))) stop('That matrix name is already used!')

print_time("Add Matrix")
getAvailableMatrices(proj)
proj <- addFeatureMatrix(input = proj, features = features_GR, matrixName = matrix_name, ceiling = ceiling_num, 
                         binarize = binarize)
getAvailableMatrices(proj)

print_time("Extract Matrix")
extract_pxc <- getMatrixFromProject(ArchRProj = proj, useMatrix = matrix_name, binarize = binarize)
if(binarize){
    saveRDS(extract_pxc,paste(sep='',outPrefix,'_sumExp_binary.rds'))
} else {
    saveRDS(extract_pxc,paste(sep='',outPrefix,'_sumExp_ceiling',ceiling_num,'.rds'))
}

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")


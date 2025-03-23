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


parser <- ArgumentParser(description='Add peaks and peaks x cells to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('peak_file', metavar='peak_file', help='peak GR file')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--extract_addToFileName', metavar='extract_addToFileName', help='if given, add to the extract filename')
parser$add_argument('--ceiling_num', metavar='ceiling_num', default=4, type='integer', help='max number of reads in pxc; ArchR default 4')
parser$add_argument('--thread_num', metavar='thread_num', default=8, type='integer', help='number of threads to use; default 8')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
parser$add_argument('--force', action='store_true', help='if given, force remake peaks; pxc is always forced')
args <- parser$parse_args()

project_dir <- args$project_dir
peak_file <- args$peak_file
extract_dir <- args$extract_dir
extract_addToFileName <- args$extract_addToFileName
ceiling_num <- args$ceiling_num
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed
force <- args$force

print_time("Argument Checking")
if(!all(file.exists(project_dir,peak_file))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(extract_dir,1)!='/') project_dir<-paste(sep="",extract_dir,'/')
if(!file.exists(extract_dir)) dir.create(extract_dir)

if(!is.null(extract_addToFileName)){
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',extract_addToFileName)
} else {
    outPrefix <- paste(sep='',extract_dir,basename(project_dir))
}
cat(paste("Using:",outPrefix,'\n'))

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_addPeakMat_args.txt")
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

print_time("Load Peaks")
peaks_GR <- readRDS(peak_file)

print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)

print_time("Add Peaks")
proj <- addPeakSet(ArchRProj = proj, peakSet = peaks_GR, force = force)

print_time("Add PxC")
getAvailableMatrices(proj)
proj <- addPeakMatrix(ArchRProj = proj, ceiling = ceiling_num, force = TRUE)
getAvailableMatrices(proj)

print_time("Extract Peaks")
extract_peaks_GR <- getPeakSet(ArchRProj = proj)
saveRDS(extract_peaks_GR,paste(sep='',outPrefix,'_peaks.rds'))

print_time("Extract PxC")
extract_pxc <- getMatrixFromProject(ArchRProj = proj, useMatrix="PeakMatrix")
saveRDS(extract_pxc,paste(sep='',outPrefix,'_pxcSumExp_ceiling',ceiling_num,'.rds'))

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")

